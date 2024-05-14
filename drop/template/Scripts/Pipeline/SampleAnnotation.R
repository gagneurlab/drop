#'---
#' title: Sample Annotation Overview
#' author:
#' wb:
#'  log:
#'   - snakemake: '`sm str(tmp_dir / "SampleAnnotation.Rds")`'
#'  params:
#'   - hpoFile: '`sm cfg.get("hpoFile")`'
#'  input: 
#'   - sampleAnnotation: '`sm sa.file`'
#'  output:
#'   - hpoOverlap: '`sm touch(cfg.getProcessedDataDir() + "/sample_anno/genes_overlapping_HPO_terms.tsv")`'
#' output:
#'   html_document:
#'    code_folding: hide
#'    code_download: TRUE
#'---

#+ echo=F
saveRDS(snakemake, snakemake@log$snakemake)

suppressPackageStartupMessages({
  library(data.table)
  library(magrittr)
  library(tidyr)
})

sa <- fread(snakemake@input$sampleAnnotation)

#'
#' Number of rows and columns in the sample annotation: `r dim(sa)`

#'
#' ## Sample annotation
DT::datatable(sa, filter = 'top')

#' ## Quality control checks

# check for duplicated rows
if(sum(duplicated(sa)) > 0){
  print("The sample annotation has the following duplicated rows. Remove them.")
  sa[duplicated(sa)]
}

#' Check for RNA_IDs without a value
if(nrow(sa[is.na(RNA_ID)]) > 0){
  print("The sample annotation has some non-existent RNA_IDs. Fill them.")
  sa[is.na(RNA_ID)]
}


#' Check for nonexistent BAM files
sa[, aux1 := file.exists(RNA_BAM_FILE)]
if(any(sa$aux1 == F)){
  print('The following BAM files do not exist: ')
  DT::datatable(sa[aux1 == F])
}

#' Check for nonexistent VCF files
if(! all(sa[,is.na(DNA_VCF_FILE)])){
  sa[, aux1 := file.exists(DNA_VCF_FILE)]
  if(any(sa$aux1 == F)){
    print('The following VCF files do not exist: ')
    DT::datatable(sa[aux1 == F])
  }
}

#' Check for RNA_IDs with more than one RNA_BAM_FILE
if(sum(duplicated(unique(sa[,.(RNA_ID, RNA_BAM_FILE)])$RNA_ID)) > 0){
  print('The following RNA_IDs and RNA_BAM_FILEs do not have a 1:1 match. Correct them.')
  DT::datatable(duplicated(unique(sa[,.(RNA_ID, RNA_BAM_FILE)])$RNA_ID))
}

#' ## Barplot with DROP groups
unique(sa[,.(RNA_ID, DROP_GROUP)])$DROP_GROUP %>% strsplit(',') %>% unlist %>%
  table %>% barplot(xlab = 'DROP groups', ylab = 'Number of samples')

# Obtain genes that overlap with HPO terms
#+ echo=F
if(!is.null(sa$HPO_TERMS) & !all(is.na(sa$HPO_TERMS)) & ! all(sa$HPO_TERMS == '')){
  sa2 <- sa[, .SD[1], by = RNA_ID]
  
  filename <- ifelse(is.null(snakemake@params$hpo_file), 
                     'https://www.cmm.in.tum.de/public/paper/drop_analysis/resource/hpo_genes.tsv.gz',
                     hpo_file)
  hpo_dt <- fread(filename)
  
  sapply(1:nrow(sa2), function(i){
    hpos <- strsplit(sa2[i, HPO_TERMS], split = ',') %>% unlist
    genes <- paste(sort(hpo_dt[HPO_id %in% hpos, hgncSymbol]), collapse = ',')
    set(sa2, i, 'HPO_matching_genes', value = genes)
  }) %>% invisible()  # don't print result
  sa2 <- sa2[, .(RNA_ID, HPO_matching_genes)]
  
  fwrite(sa2, snakemake@output$hpoOverlap,
         na = NA, sep = "\t", row.names = F, quote = F)
  }

