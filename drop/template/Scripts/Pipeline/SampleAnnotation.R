#'---
#' title: Sample Annotation Overview
#' author:
#' wb:
#'  params:
#'   - tmpdir: '`sm drop.getTmpDir()`'
#'  input: 
#'   - sampleAnnotation: '`sm config["sampleAnnotation"]`'
#' output:
#'   html_document:
#'    code_folding: hide
#'    code_download: TRUE
#'---

#+echo=F
saveRDS(snakemake, file.path(snakemake@params$tmpdir, "sample_anno_overview.snakemake"))
# snakemake <- readRDS(".drop/tmp/sample_anno_overview.snakemake")

suppressPackageStartupMessages({
  library(data.table)
  library(magrittr)
})

sa <- fread(snakemake@input$sampleAnnotation)

#' Number of rows and columns in tha sample annotation
dim(sa)
DT::datatable(sa)

# check for duplicated rows
if(sum(duplicated(sa)) > 0){
  print("The sample annotation has the following duplicated rows. Remove them.")
  sa[duplicated(sa)]
}

#' Check for NAs
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
if('DNA_VCF_FILE' %in% colnames(sa)){
  sa[, aux1 := file.exists(DNA_VCF_FILE) | is.na(DNA_ID)]
  if(any(sa$aux1 == F)){
    print('The following VCF files do not exist: ')
    DT::datatable(sa[aux1 == F])
    }
}

#' Check for RNA_IDs with more than 1 RNA_BAM_FILE
if(sum(duplicated(unique(sa[,.(RNA_ID, RNA_BAM_FILE)])$RNA_ID)) > 0){
  print('The following RNA_IDs and RNA_BAM_FILEs do not have a 1:1 match. Correct them.')
  DT::datatable(duplicated(unique(sa[,.(RNA_ID, RNA_BAM_FILE)])$RNA_ID))
}

#' Barplot with DROP groups
unique(sa[,.(RNA_ID, DROP_GROUP)])$DROP_GROUP %>% strsplit(',') %>% unlist %>%
  table %>% barplot(xlab = 'DROP groups', ylab = 'Number of samples')

