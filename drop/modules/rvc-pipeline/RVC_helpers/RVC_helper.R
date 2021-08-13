
parse_rna_dna_names <- function(l){
  rna_ind <- seq(1,length(l),3)
  dna_ind <- seq(2,length(l),3)
  file_ind <- seq(3,length(l),3)
  return(data.table("RNA_ID" = l[rna_ind],"DNA_ID" = l[dna_ind],"DNA_VCF_FILE" = l[file_ind]))
}

# helper functions for RVC overview
filter_chrm <- function(vcf_table) {
  if (grepl("chr", vcf_table$CHROM[1])){
    standard_chr <- paste0('chr', c(1:22, 'X', 'Y', 'M'))
  }
  else{
    standard_chr <- as.character(c(1:22, 'X', 'Y', 'MT'))
  }
  vcf_table <- vcf_table[CHROM %in% standard_chr]
  vcf_table[grepl("^[0-9]",CHROM),CHROM:=paste0("chr",CHROM)]
  vcf_table[grepl("^[0-9]",ID),ID:=paste0("chr",ID)]
  return(vcf_table)
}     


genotype_filter <- function(vcf_table){
  vcf_table[GT %in% c("0/1","1/0","1|0","0|1"), GT:= "0/1"]
  vcf_table[GT %in% c("1/1","1|1"), GT:= "1/1"]
  return(vcf_table[GT %in% c("0/1","1/1")])
}

extract_vcf_to_file <- function(rnaseq_vcf_file,exome_vcf_file, rna_name, wes_name){
  vcf_skip <- FALSE
  
  print(paste0('Processing ',rna_name))
  
  print(paste0("RNA VCF File: ",rnaseq_vcf_file))
  print(paste0("WES VCF File: ",exome_vcf_file))
  
  
  if(!file.exists(exome_vcf_file)){
    #WES does not exist
    print("WES VCF file does not exist. WES data will be ommitted")
    vcf_skip <- TRUE
  }
  
  if (!file.exists(rnaseq_vcf_file)) {
    paste0("RNA file does not exist. Continue")
    return(0)
  }
  
  rna_vcf <- vcfR::read.vcfR(rnaseq_vcf_file)
  if(!vcf_skip){
    wes_vcf <- vcfR::read.vcfR(exome_vcf_file)
  }
  
  # extract info from rna_vcf
  rna <- data.table(rna_vcf@fix)
  rna[,ID:=paste0(CHROM,":",POS)]
  rna[,GT:= extract.gt(rna_vcf)]
  rna[,AD:= extract.gt(rna_vcf,element = "AD")]
  rna <- separate(rna, AD, c("AD_REF", "AD_ALT"),sep = ",")
  rna[,c("QUAL","AD_REF","AD_ALT") := lapply(rna[,.(QUAL,AD_REF,AD_ALT)],as.numeric)]
  rna[, AD_ALT_FREQ := AD_ALT / (AD_ALT + AD_REF) ]
  rna[,DP:= as.numeric(extract.gt(rna_vcf,element = "DP"))]
  rna_info <- vcfR::extract_info_tidy(rna_vcf,info_fields = c("AC","AF","FS","MQ","QD"),
                                      info_types = c(AC="n",AF="n",FS="n",MQ="n",QD="n"))
  rna[,names(rna_info)[2:6]:=rna_info[,2:6]]
  # extract info from wes_vcf
  if(!vcf_skip){
    wes <- data.table(wes_vcf@fix)
    wes[,ID:= paste0(CHROM,":",POS)]
    wes[,GT:= extract.gt(wes_vcf,IDtoRowNames = F)[,wes_name]]
    wes[,DP:= extract.gt(wes_vcf,element = "DP",as.numeric = T,IDtoRowNames = F)[,wes_name]]
    wes <- wes[!duplicated(ID)]
  }
  
  
  
  # remove non-standard chromosomes
  rna <- filter_chrm(rna)
  rna <- genotype_filter(rna)
  
  if(!vcf_skip){
    wes <- filter_chrm(wes)
    wes <- genotype_filter(wes)
  }
  
  
  # Merge dataframes
  if(!vcf_skip){
    merged_df <- data.table::merge.data.table(rna,wes,by="ID",all.x=T,all.y=T, suffixes = c(".RNA",".WES"))
    print(names(merged_df))
    
    merged_df[is.na(REF.RNA),REF.RNA := REF.WES]
    merged_df <- merged_df[, .(ID,REF.RNA,ALT.RNA,ALT.WES, GT.RNA,GT.WES,QUAL.RNA,FILTER.RNA, AD_REF,AD_ALT,AD_ALT_FREQ,
                               DP.RNA,AC,AF,FS,MQ,QD,QUAL.WES,FILTER.WES,DP.WES)]
    merged_df <- merged_df %>% separate(ID, c("CHROM", "POS"),remove = F)
    merged_df <- dplyr::rename(merged_df, REF = REF.RNA)
  }
  else {
    merged_df <- rna
  }
    return(merged_df)
}

