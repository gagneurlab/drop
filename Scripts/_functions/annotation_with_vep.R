# author: christian mertes
# annotates a given vcf with vep
#

# load libraries
suppressPackageStartupMessages({
    library(VariantAnnotation)
    library(AnnotationHub)
    library(ensemblVEP)
    library(data.table)
})


SO_HASH <- list(
    "transcript_ablation" 					= c("ablation", 		1),
    "splice_donor_variant"					= c("splice_donor", 			2),
    "splice_acceptor_variant"	 			= c("splice_acceptor", 			3),
    "stop_gained" 							= c("stop", 			4),
    "frameshift_variant" 					= c("frame-shift", 		5),
    "stop_lost" 							= c("unstop", 			6),
    "start_lost"							= c("unstart",          6),
    "initiator_codon_variant"	 			= c("init", 			7),
    "transcript_amplification"				= c("amplify", 			8),
    "inframe_insertion" 					= c("ins", 				9),
    "inframe_deletion"						= c("del", 				10),
    "missense_variant"						= c("missense", 		11),
    "protein_altering_variant"				= c("missense",         11),
    "splice_region_variant"					= c("splice", 			12),
    "incomplete_terminal_codon_variant"	 	= c("incomplete", 		13),
    "stop_retained_variant" 				= c("stop_retain", 		14),
    "start_retained_variant" 				= c("start_retain", 	15),
    "synonymous_variant" 					= c("synonymous",		16),
    "coding_sequence_variant"	 			= c("coding", 			17),
    "mature_miRNA_variant" 					= c("mature_miRNA", 		18),
    "5_prime_UTR_variant" 					= c("5utr", 			19),
    "3_prime_UTR_variant" 					= c("3utr", 			20),
    "non_coding_transcript_exon_variant"	= c("ncrna-exon", 		21),
    "intron_variant" 						= c("intron", 			22),
    "NMD_transcript_variant"		 		= c("NMD", 				23),
    "non_coding_transcript_variant"			= c("ncrna-transcript", 24),
    "downstream_gene_variant"		 		= c("downstream", 		25),
    "feature_elongation" 					= c("feature_e", 		26),
    "feature_truncation" 					= c("feature_t", 		27),
    "intergenic_variant" 					= c("intergenic", 		28),
    "regulatory_region_variant"				= c("regulation", 		29),
    "TF_binding_site_variant" 				= c("binding", 			30),
    "upstream_gene_variant" 				= c("upstream", 		31)
)

#'
#' sets the needed parameters to run the VEP tool from ENSEMBL 
#' 
get_vep_params <- function(version=max(unlist(currentVEP())), num_forks=4, 
                           local=FALSE, vcfFile=tempfile("vep_anno-", fileext = ".vcf.gz")){
    
    # define vep pathes
    addVEP2Path(version)
    module_path   <- "/opt/modules"
    vep_module_id <- file.path("i12g/ensembl-vep", version)
    vep_cache_dir <- file.path(module_path, vep_module_id, "cachedir")
    
    # set environment variables
    
    
    # create new vep param object
    vep_param <- VEPFlags(version)
    flags(vep_param)$output_file     <- vcfFile
    scriptPath(vep_param)            <- system("which vep", intern=TRUE)
    flags(vep_param)$v               <- TRUE
    flags(vep_param)$vcf             <- TRUE
    flags(vep_param)$compress_output <- "bgzip"
    flags(vep_param)$minimal         <- TRUE
    flags(vep_param)$allele_number   <- TRUE
    flags(vep_param)$no_stats        <- FALSE
    flags(vep_param)$stats_file      <- paste0(vcfFile, "_summary.html")
    flags(vep_param)$force_overwrite <- TRUE
    
    # set basic vep options
    flags(vep_param)$verbose    <- TRUE
    flags(vep_param)$everything <- TRUE
    flags(vep_param)$fork       <- num_forks
    flags(vep_param)$assembly   <- "GRCh37"
    
    # vep data access options
    flags(vep_param)$db_version <- version
    flags(vep_param)$merged     <- TRUE
    flags(vep_param)$user       <- "anonymous"
    #flags(vep_param)$port       <- 5306 # default port
    flags(vep_param)$port       <- 3337 # used for GRCh37
    flags(vep_param)$host       <- "ensembldb.ensembl.org"
    
    # vep cache options
    if(flags(vep_param)$merged){
        flags(vep_param)$database <- FALSE
    }
    flags(vep_param)$cache       <- TRUE
    flags(vep_param)$offline     <- FALSE
    flags(vep_param)$dir         <- vep_cache_dir
    flags(vep_param)$dir_cache   <- vep_cache_dir
    flags(vep_param)$dir_plugins <- file.path(vep_cache_dir, "Plugins")
    
    # vep advanced options
    flags(vep_param)$buffer_size <- 10000
    
    # output options
    flags(vep_param)$sift         <- 's'
    flags(vep_param)$polyphen     <- 's'
    flags(vep_param)$total_length <- TRUE
    flags(vep_param)$numbers      <- TRUE
    
    # vep annotation options
    flags(vep_param)$symbol      <- TRUE
    flags(vep_param)$hgvs        <- TRUE
    flags(vep_param)$ccds        <- TRUE
    flags(vep_param)$uniprot     <- TRUE
    flags(vep_param)$xref_refseq <- TRUE
    flags(vep_param)$af          <- TRUE
    flags(vep_param)$max_af      <- TRUE
    flags(vep_param)$af_exac     <- TRUE
    flags(vep_param)$af_gnomad   <- TRUE
    flags(vep_param)$pubmed      <- TRUE
    flags(vep_param)$canonical   <- TRUE
    flags(vep_param)$biotype     <- TRUE
    flags(vep_param)$failed      <- 1
    
    # add CADD
    flags(vep_param)$plugin <- "CADD,/s/genomes/human/hg19/CADD/v1.3/whole_genome_SNVs.tsv.gz,/s/genomes/human/hg19/CADD/v1.3/InDels.tsv.gz" # --plugin MMSplice"
    
    #flags(vep_param)$terms       <- "ENSEMBL"
    
    # use only default chromosomes for annotation (leave the unaligned once blank
    flags(vep_param)$chr <- paste0("chr", c(1:22, "X","Y","M"), collapse=",")
    
    # return vep param object
    return(vep_param)
}

#'
#' calculates the mtype (snp, ins, del, strucVar) from the given data
#' 
get_mtype_annotation <- function(data){
    
    tmp_dt <- data.table(alt=as.character(data[,alt]), ref=as.character(data[,ref]))
    tmp_dt[,refn:=nchar(ref)]
    tmp_dt[,altn:=nchar(alt)]
    
    dc_mtype <- tmp_dt[,list(mtype=
                                 ifelse(refn == 1 & altn == 1, "snp",
                                        ifelse(refn == 1 & altn >  1, "ins", 
                                               ifelse(refn >  1 & altn == 1, "del", 
                                                      ifelse(refn >  1 & altn >  1, "strucVar",
                                                             paste("should not happen! in get_mtype_annotation '", ref, ":", alt, "'"))
                                               ))))]
    
    return(as.factor(dc_mtype[,mtype]))
}

#'
#' returns the given vector if the length of the given object is not the same as the check_size.
#' Else it will return NA 
#' 
get_vector_or_na <- function(vector, check_size = 0){
    if(length(vector) == check_size){
        return(NA)
    }
    return(vector)
}

#'
#' creates the vcf data table from the given vcf object
#' 
get_vcf_data_table <- function(sample_id, vcf_obj){
    
    dc_chr   <- gsub("^chr", "", seqnames(vcf_obj))
    
    if(!any(width(ref(vcf_obj)) != 0)){
        stop("What happened here! Why do I have this here? For indels?")
        dc_ref <- end(ranges(vcf_obj))
        dc_alt <- NA
    } else {
        dc_ref <- as.factor(as.data.frame(ref(vcf_obj))$x)
        # take the first variant if we have two
        dc_alt <- data.table(as.data.frame(alt(vcf_obj)))
        dc_alt <- as.factor(dc_alt[!duplicated(group),value])
    }
    dc_mtype <- get_mtype_annotation(data.table(ref=dc_ref,alt=dc_alt))
    
    if(is.null(info(vcf_obj)$DP4)){
        dc_dp41 <- dc_dp42 <- dc_dp43 <- dc_dp44 <- NULL
    } else {
        dp4_matrix <- t(matrix(data.table(as.data.frame(info(vcf_obj)$DP4))[,value],nrow=4))
        dc_dp41  <- dp4_matrix[,1]
        dc_dp42  <- dp4_matrix[,2]
        dc_dp43  <- dp4_matrix[,3]
        dc_dp44  <- dp4_matrix[,4]
    }
    
    # formate pl column
    pl_dt <- data.table(
        pl1=sapply(geno(vcf_obj)$PL, function(x) x[1]), 
        pl2=sapply(geno(vcf_obj)$PL, function(x) x[2]), 
        pl3=sapply(geno(vcf_obj)$PL, function(x) x[3])
    )
    dc_pl <- pl_dt[,paste("[",pl1,"|",pl2,"|",pl3,"]", sep="")]
    
    stable <- data.table(
        # sample information
        sample  = as.factor(sample_id),
        
        # variant definition
        var_id  = as.factor(1:length(rownames(vcf_obj))),
        chr     = as.factor(dc_chr),
        pos     = as.integer(start(vcf_obj)),
        ref     = as.factor(dc_ref),
        alt     = as.factor(dc_alt),
        gt      = as.factor(geno(vcf_obj)$GT),
        gq      = as.integer(geno(vcf_obj)$GQ),
        
        # quality annotation
        qual    = as.integer(qual(vcf_obj)),
        dp      = as.integer(info(vcf_obj)$DP),
        # for compatibilty with deletions from Dindel
        dp41    = as.integer(get_vector_or_na(dc_dp41)),
        dp42    = as.integer(get_vector_or_na(dc_dp42)),
        dp43    = as.integer(get_vector_or_na(dc_dp43)),
        dp44    = as.integer(get_vector_or_na(dc_dp44)),
        mq      = as.integer(get_vector_or_na(info(vcf_obj)$MQ)),
        sp      = as.integer(get_vector_or_na(geno(vcf_obj)$SP)),
        pl      = as.factor(get_vector_or_na(dc_pl, 2)),
        mtype   = as.factor(dc_mtype)
    )
    
    if(length(vcf_obj) != nrow(stable)){
        stop("Error in extraction of rows for data table!")
    }
    
    return(stable)
}


simplify_so_terms <- function(so_terms, sample){
    
    # simplify SO terms and add a severity index to it
    FUN <- function(x){
        if(is.na(x)){
            return(c(NA, NA))
        }
        # retrieve simplified term from SO_HASH
        so_hash_entry <- unlist(SO_HASH[strsplit(x, "&")[[1]]])
        if (is.null(so_hash_entry)) {
            msg <- paste0("sample ", sample, ": term ", x, " not found when simplifying SO terms")
            message(msg)
            write(msg, '~/Downloads/vep2dt.errors', append = T)
            return(c(x,x)) # keep unsimplified term
        }
        
        matrix_so <- t(matrix(so_hash_entry,nrow=2))
        if(dim(matrix_so)[1] == 1){	
            return(as.vector(t(matrix_so)))
        } else {
            return(matrix_so[which.min(matrix_so[,2]),])
        }
    }
    so_terms <- as.factor(so_terms)
    results <- t(matrix(sapply(levels(so_terms), FUN),nrow=2))
    return(list(
        dc_mstype=as.factor(results[as.numeric(so_terms),1]),
        dc_severe=as.integer(results[as.numeric(so_terms),2])
    ))
}

#'
#' all annotated frequencies by default
#' 
get_frequencies_from_vep <- function(vep_obj, db = c('ExAC', 'gnomAD')){
    
    #' description of MAF on page 16 https://m.ensembl.org/info/docs/tools/vep/online/VEP_web_documentation.pdf
    #' *AF: 1000 Genomes or ESP
    #' gnomad_*_AF: gnomAD
    #' ExAC_*_AF: ExAC
    #' MAX_AF: Maximum observed allele frequency in 1000 Genomes, ESP and gnomAD
    #' MAX_AF_POPS: Populations in which maximum allele frequency was observed
    
    maf_cols <- grep('AF', colnames(mcols(vep_obj)), value = T)
    maf_cols <- grep('MAX_AF_POPS', maf_cols, value = T, invert = T)
    for (d in c('ExAC', 'gnomAD')) {
        if (!d %in% db)
            maf_cols <- grep(d, maf_cols, value = T, invert = T)
    }
    maf_dt <- as.data.table(mcols(vep_obj)[, maf_cols])
    maf_dt <- maf_dt[, lapply(.SD, as.double)]
    
    # AF         <- as.double(vep_obj$AF)
    # gnomAD_AF <- as.double(vep_obj$gnomAD_AF)
    # gnomAD_NFE_AF    <- as.double(vep_obj$gnomAD_NFE_AF)
    # gnomAD_AFR_AF     <- as.double(vep_obj$gnomAD_AFR_AF)
    # MAX_AF    <- as.double(vep_obj$MAX_AF)
    
    maf_dt
}

get_vep_annotation_data_table <- function(vep_obj, sample = 'sample'){
    
    # genename annotation
    dc_hgncid  <- ifelse(as.character(vep_obj$SYMBOL_SOURCE) == "HGNC", as.character(vep_obj$SYMBOL), NA)
    dc_uniprot <- gsub("_HUMAN$", "", as.character(vep_obj$SWISSPROT))
    dc_ensgid  <- gsub("ENSG0*", "", vep_obj$Gene)
    
    # transcipt and exon annotation
    dc_enstid  <- gsub("ENST0*", "", vep_obj$Feature)
    dc_enstid[grepl("^(ENSR|MA|PB|(NC|NM|NR|XM|XR)_)", dc_enstid, perl = TRUE)] <- NA 
    dc_tunum   <- gsub("(-[0-9]+)*/[0-9]*$", "", vep_obj$EXON)
    
    # consequences and severeness 
    dc_consequences_list <- simplify_so_terms(vep_obj$Consequence, sample)
    dc_mstype  <- dc_consequences_list[["dc_mstype"]]
    dc_severe   <- dc_consequences_list[["dc_severe"]]
    
    # mafs
    maf_table <- get_frequencies_from_vep(vep_obj, db = "gnomAD")
    
    # mutation references like pubmed paper or SNPdb aka rsid
    dc_rsid    <- gsub("^(rs\\d+)([^\\d].*)?$", "\\1", as.character(vep_obj$Existing_variation))
    dc_rsid    <- ifelse(grepl("^rs[0-9]*$", dc_rsid), dc_rsid, NA)
    
    # create data.table
    vep_table <- data.table(
        
        # variation id
        var_id   = as.factor(vep_obj$VCFRowID),
        
        # genename annotation
        hgncid   = as.factor(dc_hgncid),
        sift1    = as.double(gsub(".*\\(|\\)", "", vep_obj$SIFT)),
        pph1     = as.double(gsub(".*\\(|\\)", "", vep_obj$PolyPhen)),
        uniprot  = as.factor( dc_uniprot),
        ucsc     = NA,
        ccds     = as.factor( vep_obj$CCDS),
        refseq   = as.factor( vep_obj$RefSeq),
        noccds   = as.logical(is.na(vep_obj$CCDS)),
        norefseq = as.logical(is.na(vep_obj$RefSeq)),
        
        # vep annotation
        ensgid   = as.integer(dc_ensgid),
        enstid   = as.integer(dc_enstid),
        exonic   = as.logical(!is.na(vep_obj$EXON)),
        tunum    = as.integer(dc_tunum),
        sdistl   = as.integer(vep_obj$DISTANCE),
        sdistr   = as.integer(vep_obj$DISTANCE),
        mstype   = as.factor( dc_mstype),
        # af       = as.double( maf_table[,af]),
        # gnomad_maf = as.double(maf_table[,gnomad_maf]),
        # nfe_maf    = as.double(maf_table[,nfe_maf]),
        # aa_maf     = as.double(maf_table[,aa_maf]),
        # max_maf    = as.double(maf_table[,max_maf]),
        rsid     = as.factor( dc_rsid),
        rspm     = NA,
        rsg5     = NA,
        pubmed   = as.factor( vep_obj$PUBMED),
        severe   = as.integer(dc_severe),
        
        CADD_phred = vep_obj$CADD_PHRED,
        CADD_raw = vep_obj$CADD_RAW
    )
    vep_table <- cbind(vep_table, maf_table)
    
    return(vep_table)
}

parse_patient_id <- function(file_name){
    file_name <- as.character(file_name)
    sample_id <-
        # helmholtz data structure
        ifelse(grepl("ontarget.varfilter.dbSNP.plus.checked.vcf$", file_name, perl = T ),
               basename(dirname(dirname(dirname(file_name)))),
               # klein data structure
               ifelse(grepl("((scn)|(ibd)|(sy))\\d+[a-w]{2}", basename(file_name), perl = T), 
                      gsub("[_.].*", "", basename(file_name), perl = T),
                      # no id extracted
                      paste("can't identify sample id from file name: '", file_name, "'")
               ))
    
    return(sample_id)
}

#' 
#' annotate vcf and save to file
#' 
run_vep_annotation <- function(vcf_file, vepVcfFile, num_forks, vep_version=DEFAULT_VEP_VERSION){
    
    #vepVcfFile <- paste0(outDir, '/', sample, '-', file_idx, '-annotation_vep.vcf.gz')
    
    message("Annotating file ", vcf_file, " to ", vepVcfFile, print_stack=FALSE)
    
    vep_param <- get_vep_params(vep_version, num_forks=num_forks, vcfFile=vepVcfFile)
    ensemblVEP(vcf_file, vep_param)  # The vep_param already contains the output file
}

#' TODO: sample redundant
#' combine vcf and vep tables
#' 
combine_vcf_vep <- function(sample, vcf_obj, vep_obj, minQUAL=20, num_forks) {
    
    # create data.tables for sample
    message("create vcf data table ...")
    vcf_data_table <- get_vcf_data_table(sample, vcf_obj)
    message("create vep data table ...")
    vep_data_table <- get_vep_annotation_data_table(vep_obj, sample)
    
    # combine both data.tables
    setkey(vcf_data_table, var_id)
    setkey(vep_data_table, var_id)
    annotated_data_table <- vcf_data_table[vep_data_table]
    
    # filter by severity and consensus
    setkey(annotated_data_table, var_id, noccds, severe, norefseq)
    uniq_annotated_data_table <- annotated_data_table[!duplicated(annotated_data_table[,var_id])]
    
    if(dim(vcf_obj)[1] != nrow(uniq_annotated_data_table)){
        message("lost some variations during mapping for sample '", sample, "'!!!")
    }
    return(uniq_annotated_data_table)
    
}

run_annotation_and_save_it <- function(sample, vcf_file, outDir, num_forks,
                                       overrideAll=FALSE, overrideRdata=FALSE, 
                                       minQUAL=20, 
                                       vep_version=DEFAULT_VEP_VERSION,
                                       rdsFile=paste0(outDir, '/', sample, '_annotated_data_table.rds')){
    
    if(grepl("^can't identify sample id", sample)){
        message("\n", sample)
        return()
    }
    
    message("Working on patient ", sample, " ... ")
    
    # if more files exist combine the output and save it to one file
    combined_annotated_data_table <- NULL
    
    for(file in vcf_files){
        file_idx <- which(file == vcf_files)
        
        # skip this file if it has no variants
        if(nrow(readVcf(file, genome = 'hg19')) < 1){
            message("File '", file, "' has no variants! It will be skipped!")
            next
        }
        
        message("Run for file ", file, " with index ", file_idx, ' ... ',
                  print_stack=FALSE)
        
        # files which will  be created
        vepVcfFile <- paste0(outDir, '/', sample, '-', file_idx, 
                             '-annotation_vep.vcf.gz')
        
        if( overrideAll || ! file.exists(vepVcfFile) ){
            # annotate file
            message("Read vcf & annotate it with VEP ... ", print_stack=FALSE)
            vep_param <- get_vep_params(vep_version, num_forks=num_forks, 
                                        vcfFile=vepVcfFile)
            ensemblVEP(file, vep_param)  # The vep_param already contains the output file
            
        } else {
            message("using existing annotation file ...")
        }
        
        if( overrideAll || overrideRdata || ! file.exists(rdsFile) ) {
            vcf_obj <- readVcf(vepVcfFile, "hg19")
            # discard any mutations with less minQUAL
            vcf_obj <- vcf_obj[rowRanges(vcf_obj)$QUAL >= minQUAL]
            vep_obj <- parseCSQToGRanges(vcf_obj, VCFRowID = rownames(vcf_obj))
            
            # create data.tables for sample
            message("create vcf data table ...")
            vcf_data_table <- get_vcf_data_table(sample, vcf_obj)
            message("create vep data table ...")
            vep_data_table <- get_vep_annotation_data_table(vep_obj) #, num_forks)
            
            # combine both data.tables
            setkey(vcf_data_table, var_id)
            setkey(vep_data_table, var_id)
            annotated_data_table <- vcf_data_table[vep_data_table]
            
            # take only most severe annotation
            setkey(annotated_data_table, var_id, noccds, severe, norefseq)
            uniq_annotated_data_table <- annotated_data_table[!duplicated(annotated_data_table[,var_id])]
            
            if(dim(vcf_obj)[1] != nrow(uniq_annotated_data_table)){
                message("lost some variations during mapping for sample '", sample, "'!!!")
            }
            
            # combine multiple files into one data.table
            if(is.null(combined_annotated_data_table)){
                combined_annotated_data_table <- uniq_annotated_data_table
            } else {
                combined_annotated_data_table <- rbind(combined_annotated_data_table, uniq_annotated_data_table)
            }
        }
    }
    
    # save object
    message('save r-object for patient ', sample, " ... ")
    saveRDS(combined_annotated_data_table, rdsFile)
    return(rdsFile)
}

#'
#' annotates the given cnv file or data table with ensembl's vep tool
#' 
run_vep_annotation_for_cnv_calls <- function(cnv_calls, output_dir, num_forks = 4){
    # create dir if not there yet
    if(!file.exists(output_dir)){
        dir.create(output_dir, TRUE, TRUE)
    }
    
    # check if cnv calls is provided as a data table
    if(any(class(cnv_calls) %in% "data.table")){
        vep_struc_var_in_file <- file.path(output_dir, "cnv_calls_vep_in.vcf")
        vep_vcf_out_file <- paste0(output_dir, "/cnv_call_vep_out.vcf.gz")
        rds_cnv_file <- paste0(output_dir, "/annotation_cnv_call.rds") 
        
        # write vep input file
        write_vcf_like_cnv_file(cnv_calls, vep_struc_var_in_file)
        bed_like_obj <- cnv_calls
        cnv_calls1 <- "given cnv data table"
        
        # check if a vcf is provided or a bed like file
    } else if(grepl("*.vcf(.gz)?$", cnv_calls)){
        base_name_file <- gsub(".vcf(.gz)?$", "", basename(cnv_calls), perl = TRUE)
        vep_struc_var_in_file <- file.path(output_dir, paste0(base_name_file, ".vcf"))
        vep_vcf_out_file <- paste0(output_dir, "/vep_vcf_out_", base_name_file, ".vcf.gz")
        rds_cnv_file <- paste0(output_dir, "/annotation_", base_name_file, ".rds") 
        
        # write vep input file and remove deletions with length -1 
        # TODO (error in vep ==> not equal number of lines after annotations!)
        if(grepl(".vcf.gz", cnv_calls)){
            grep_cmd <- "zgrep"
        } else {
            grep_cmd <- "grep"
        }
        message("Removed deletions with length -1 from cnv lumpy calls from file '", cnv_calls, "'.\n",
                  paste(system(paste(grep_cmd, " '^chr.*SVTYPE=DEL;.*SVLEN=-1;'", cnv_calls), intern = TRUE), collapse = "\n"),
                  print_stack = FALSE
        )
        system(paste(grep_cmd, "-v '^chr.*SVTYPE=DEL;.*SVLEN=-1;'", cnv_calls, " > ", vep_struc_var_in_file))
        
        # read the bed like file
        message("read bedlike cnv file '", base_name_file, "' ...")
        bed_like_obj <- read_vcflike_cnv_file(vep_struc_var_in_file, TRUE)
        cnv_calls1 <- cnv_calls
        
    } else { 
        # needed file names 
        base_name_file <- basename(cnv_calls)
        vep_struc_var_in_file <- paste0(output_dir, "/vep_in_", base_name_file, ".bed")
        vep_vcf_out_file <- paste0(output_dir, "/vep_vcf_out_", base_name_file, ".vcf.gz")
        rds_cnv_file <- paste0(output_dir, "/annotation_", base_name_file, ".rds") 
        
        # read the bed like file
        message("read bedlike cnv file ", base_name_file)
        bed_like_obj <- read_bedlike_cnv_file(cnv_calls, TRUE)
        
        # write vep input file
        write_vcf_like_cnv_file(bed_like_obj, vep_struc_var_in_file)
        cnv_calls1 <- cnv_calls
    }
    
    # run vep
    message("run vep from ensemble with file ", vep_struc_var_in_file)
    vep_param <- get_vep_params(80, num_forks, local = FALSE)
    ensemblVEP(vep_struc_var_in_file, vep_param)
    system(paste("cat ",	flags(vep_param)$output_file, " | gzip > ", vep_vcf_out_file, sep=""), wait = TRUE)
    
    # read vcf file 
    message("read vcf file from vep output file ", vep_vcf_out_file)
    vcf_obj <- readVcf(vep_vcf_out_file, "hg19")
    
    # add geno type informations
    vcf_obj <- add_genotype_informations(vcf_obj, bed_like_obj$genotype)
    
    # get vcf data table
    vcf_data_table <- get_vcf_data_table("bed_like_sample", vcf_obj, num_forks = num_forks)
    vcf_data_table <- add_bed_like_informations(vcf_data_table, bed_like_obj)
    
    # add dgv frequencies to it
    vcf_data_table <- add_dgv_informations(vcf_obj, vcf_data_table, bed_like_obj, num_forks)
    
    # get vep data table
    vep_obj <- parseCSQToGRanges(vcf_obj, VCFRowID = rownames(vcf_obj))
    vep_data_table <- get_vep_annotation_data_table(vep_obj, num_forks)
    
    # combine both data.tables
    message("combine vep with vcf file informations ...")
    setkey(vcf_data_table, var_id)
    setkey(vep_data_table, var_id)
    annotated_data_table <- vcf_data_table[vep_data_table]
    
    # put hgncid into alt column
    annotated_data_table[,alt:=hgncid]
    
    # take only most severe annotation
    setkey(annotated_data_table, var_id, alt, noccds, severe, norefseq)
    uniq_annotated_data_table <- annotated_data_table[!duplicated(annotated_data_table[,list(var_id,hgncid)])]
    
    # save rds object
    message('save r-object for cnv file "', cnv_calls1 ,'" to "', rds_cnv_file, '" ... ')
    saveRDS(uniq_annotated_data_table, rds_cnv_file)
    
    # remove tmp files
    file.remove(flags(vep_param)$output_file)
    file.remove(vep_struc_var_in_file)
    
    # return the rds file name
    return(rds_cnv_file)
}

#'
#' 
add_dgv_informations <- function(vcf_obj, vcf_data_table, bed_like_obj, num_forks = 4,
                                 dgv_file = "/s/genomes/human/hg19/dgv/GRCh37_hg19_variants_2014-10-16_loss_n_deletions.txt"){
    # read dgv data and convert the counts into integers
    dgv_data <- read_bedlike_cnv_file(dgv_file)
    mcols(dgv_data)$pop_size <- as.integer(as.character(mcols(dgv_data)$pop_size))
    mcols(dgv_data)$loss_counts <- as.integer(as.character(mcols(dgv_data)$loss_counts))
    mcols(dgv_data)$gain_counts <- as.integer(as.character(mcols(dgv_data)$gain_counts))
    
    # correct the start/end postion of the rowData
    range_vcf <- rowRanges(vcf_obj)
    end(range_vcf) <- bed_like_obj[,end]
    
    # find overlaps
    overlaps <- findOverlaps(range_vcf, dgv_data)
    
    # get only overlaps with at least 60 percent overlap
    overlap_width <- width(ranges(overlaps, ranges(range_vcf), ranges(dgv_data)))
    good_overlaps <- overlap_width / width(range_vcf[overlaps@queryHits]) > 0.6
    overlaps <- overlaps[good_overlaps]
    
    FUN <- function(q_id, ov, dgv_d){
        dgv_mcols <- mcols(dgv_d[ov@subjectHits[ov@queryHits == q_id]])
        pop_size  <- sum(dgv_mcols$pop_size)
        del_size  <- sum(dgv_mcols$loss_counts)
        gain_size <- sum(dgv_mcols$gain_counts)
        return(c(pop_size = pop_size, loss_counts = del_size, gain_counts = gain_size))
    }
    
    dgv_freq_results <- mclapply(1:dim(vcf_obj)[1], mc.cores = num_forks,
                                 FUN = FUN, ov = overlaps, dgv_d = dgv_data
    )
    
    # add frequencies to the table
    vcf_data_table[,dgv_pop_size:=sapply(dgv_freq_results, function(x) x['pop_size'])]
    vcf_data_table[,dgv_gain_counts:=sapply(dgv_freq_results, function(x) x['gain_counts'])]
    vcf_data_table[,dgv_loss_counts:=sapply(dgv_freq_results, function(x) x['loss_counts'])]
    
    return(vcf_data_table)
}

#'
#' 
add_bed_like_informations <- function(vcf_data_table, bed_like_obj){
    if(any(vcf_data_table[,chr] != gsub("^chr", "", bed_like_obj[,chr]))){
        stop("vcf file differeces from bed like object in the chromosome column!")
    }
    if(any(vcf_data_table[,pos] != bed_like_obj[,start])){
        stop("vcf file differeces from bed like object in the position/start column!")
    }
    
    vcf_data_table[,sample:=bed_like_obj[,sample_id]]
    vcf_data_table[,ref:=bed_like_obj[,end]]
    
    return(vcf_data_table)
}

#'
#' 
add_genotype_informations <- function(vcf_obj, genotype = rep(NA, dim(vcf_obj)[1]), 
                                      geno_quality = rep(NA, dim(vcf_obj)[1]), geno_scores = rep(NA, dim(vcf_obj)[1])){
    
    # add sample informations
    colData(vcf_obj) <- DataFrame(row.names = "sample_tmp", Samples = "sample_tmp")
    
    # create geno type matrices
    if(is.vector(genotype) || is.factor(genotype)){
        genotype <- matrix(genotype, ncol = 1)
    }
    if(is.vector(geno_quality) || is.factor(geno_quality)){
        geno_quality <- matrix(geno_quality, ncol = 1)
    }
    if(is.vector(geno_scores) || is.factor(geno_scores)){
        geno_scores <- matrix(geno_scores, ncol = 1)
    }
    if(any(genotype == "HET" || genotype == "HOM" || genotype == "DUP")){
        genotype <- ifelse(genotype == "HOM", "0/1", 
                           ifelse(genotype == "HET", "1/1",
                                  ifelse(genotype == "DUP", "./.", NA)))
    }
    
    # add header informations
    geno(header(vcf_obj)) <- DataFrame(row.names = c("GT", "GQ", "SP"),
                                       Number = c('.', '.', '.'), 
                                       Type = c('String', 'String', 'String'), 
                                       Description = c('genotype', 'genotype quality', 'SP column')
    )
    
    # add it to the vcf object  
    geno(vcf_obj) <- SimpleList(GT = genotype, GQ = geno_quality, SP = geno_scores)
    
    # add additional variables	
    vcf_length <- dim(vcf_obj)[1]
    if(!any(row.names(info(header(vcf_obj))) %in% "DP")){
        info(header(vcf_obj)) <- rbind(info(header(vcf_obj)),
                                       DataFrame(row.names = "DP", Number = '.', Type = 'String', Description = "Depth"))
        info(vcf_obj)$DP <- rep(NA, vcf_length)
    }
    if(!any(row.names(info(header(vcf_obj))) %in% "MQ")){
        info(header(vcf_obj)) <- rbind(info(header(vcf_obj)),
                                       DataFrame(row.names = "MQ", Number = '.', Type = 'String', Description = "Mapquality"))
        info(vcf_obj)$MQ <- rep(NA, vcf_length)
    }
    
    return(vcf_obj)
}

addVEP2Path <- function(version, perlVer="5.24.0", samtoolVer="1.3.1"){
    addSysEnv(paste0("/opt/modules/i12g/ensembl-vep/", version, "/bin"))
    addSysEnv(paste0("/opt/modules/i12g/ensembl-vep/", version, "/modules/htslib"))
    addSysEnv(paste0("/opt/modules/i12g/ensembl-vep/", version, "/modules"), "PERL5LIB")
    addSysEnv(paste0("/opt/modules/i12g/perl/", perlVer, "/bin"))
    addSysEnv(paste0("/opt/modules/i12g/perl/", version, "/perl5lib/lib/perl5"), "PERL5LIB")
    addSysEnv(paste0("/opt/modules/i12g/samtools/", samtoolVer, "/bin"))
}

addSysEnv <- function(value, var="PATH", sep=":"){
    env <- Sys.getenv(var)
    if(!grepl(value, env)){
        newenv <- paste(value, env, sep=sep)
        cmd <- paste0("Sys.setenv(", var, "=newenv)")
        eval(parse(text=cmd))
    }
}

debug_test_VEP_command <- function(){
    params <- get_vep_params(version=94)
    ensemblVEP("/s/project/mitoMultiOmics/raw_data/helmholtz/33254/exomicout/paired-endout/stdFilenames/33254.vcf.gz", param=params, verbose=TRUE)    
}