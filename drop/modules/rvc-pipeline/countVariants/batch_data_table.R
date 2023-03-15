#'---
#' title: "RNA Variant Calling Data Table"
#' author: nickhsmith
#' wb:
#'  log:
#'    - snakemake: '`sm str(tmp_dir / "RVC" / "{dataset}" / "{annotation}_RVC_data_table.Rds")`'
#'  input:
#'   - configParams: '`sm os.path.join(
#'                        cfg.processedDataDir,
#'                        "rnaVariantCalling/params/config/rnaVariantCalling_config.tsv")`'
#'   - annotatedVCF: '`sm os.path.join(
#'                        cfg.processedResultsDir,
#'                        "rnaVariantCalling/batch_vcfs", "{dataset}",
#'                        "{dataset}_{annotation}.annotated.vcf.gz")`'
#'  output:
#'   - data_table: '`sm os.path.join(
#'                        cfg.processedResultsDir,
#'                        "rnaVariantCalling/data_tables", "{dataset}",
#'                        "{dataset}_{annotation}_data_table.Rds")`'
#'  type: script
#'---

#+ echo=FALSE
library(data.table)
library(VariantAnnotation)
library(tMAE)
library(dplyr)
library(GenomicScores)

saveRDS(snakemake, snakemake@log$snakemake)
####
# Helper functions
####
.clean_vars <- function(x){
  x[which(x == "1/1" | x == "1|1")] <- "1/1"
  x[which(x == "1/0" | x == "1|0" | x == "0/1" | x == "0|1")] <- "0/1"
  x[which(x == "./." | x == ".|." | x == "0/0 "| x == "0|0")] <- "0/0"
  x[grepl("QD|FS|SnpCluster",x)] <- "Seq_filter"

  return(x)
}

# extract gtf_label by using a regex to match the target label and then take the first field after that
.extract_info <- function(vcf,gtf_label){
  gene_info <- gsub('"', "", info(vcf)$GENE)
  target <- sapply(gene_info,function(x){
    backhalf <- gsub(paste0(".*",gtf_label," (.+)"), "\\1", x)
    target <- strsplit(backhalf," ")[[1]][1]
  })
}

vcffile <- open(VcfFile(snakemake@input$annotatedVCF,yieldSize = snakemake@config$rnaVariantCalling$yieldSize)) # read the batch vcf
res_final <- data.table()

#while there are rows to read in. Process the vcf file until nrow is 0.
while(nrow(vcf <- readVcf(vcffile,param =  ScanVcfParam(fixed="FILTER",geno="GT")))) {
  canonical_chr <- c(paste0("chr",c(1:22,"X","Y","M")),1:22,"X","Y","MT")
  vcf <- vcf[seqnames(vcf) %in% canonical_chr]

  vcf_dt <- as.data.table(geno(vcf)$GT)
  vcf_dt[,FILTER := vcf@fixed$FILTER]
  vcf_dt[,VARIANT := names(vcf)]
  vcf_dt[,GENE_ID := .extract_info(vcf,"gene_id")]
  vcf_dt[,GENE_NAME := .extract_info(vcf,"gene_name")]

  dt <- as.data.table(lapply(vcf_dt, .clean_vars))

  # calculate variant frequency within the cohort
  max_var_freq_cutoff <- snakemake@config$rnaVariantCalling$maxVarFreqCohort
  dt$cohortFreq <- apply(dt[,1:ncol(vcf)],1,function(x){
    1-sum(x == "0/0")/length(x) })


  if (snakemake@config$rnaVariantCalling$addAF){
    #tMAE code for adding gnomad frequency
    max_af_cutoff <- snakemake@config$rnaVariantCalling$maxAF
    pops <- c('AF', 'AF_afr', 'AF_amr', 'AF_eas', 'AF_nfe', 'AF_popmax')
    maf_dt <- add_gnomAD_AF(granges(vcf),genome_assembly = snakemake@config$genomeAssembly) %>% as.data.table()

    # Compute the MAX_AF based on all provided population columns
    # return -1 if only NAs are present (to avoid a warning)
    maf_dt$MAX_AF <- apply(maf_dt[, ..pops], 1, 
                      FUN=function(x){ max(x, -1, na.rm=TRUE) })

    # Replace Inf/-1 with NA
    maf_dt[is.infinite(MAX_AF) | MAX_AF == -1, MAX_AF := NA]

    res <- cbind(dt,maf_dt[,"MAX_AF"])

    # Label variants based on their MAX_AF
    res[FILTER == "PASS",FILTER:= "PASS_common"]
    res[MAX_AF <= max_af_cutoff & cohortFreq <= max_var_freq_cutoff & FILTER == "PASS_common",FILTER:= "PASS_rare"]

    #Reorder FILTER
    res$FILTER <- factor(res$FILTER,levels = c("Seq_filter","Mask","minALT","Mask;minALT","PASS_common","PASS_rare"))

  } else{
    res <- dt
    #Reorder FILTER
    res$FILTER <- factor(res$FILTER,levels = c("Seq_filter","Mask","minALT","Mask;minALT","PASS"))
    res$MAX_AF <- NA
  }
  res_final <- rbind(res_final,res)
}
 
close(vcffile) 
setcolorder(res_final,c("VARIANT","GENE_ID","GENE_NAME","FILTER","MAX_AF","cohortFreq"))
saveRDS(res_final,snakemake@output$data_table) # save res to a data_table object
