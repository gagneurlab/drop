#'---
#' title: Sample Annotation Overview
#' author:
#' wb:
#'  params:
#'   - tmpdir: '`sm drop.getTmpDir()`'
#'  input: 
#'   - sampleAnnotation: '`sm config["sampleAnnotation"]`'
#'  output:
#'   - done: '`sm parser.getProcDataDir() + "/sample_anno/sample_anno.done"`'
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

# Obtain genes that overlap with HPO terms
#+echo=F
if(!is.null(sa$HPO_TERMS)){
  sa2 <- sa[, .SD[1], by = RNA_ID]
  hpo_dt <- fread('https://i12g-gagneurweb.in.tum.de/public/paper/drop_analysis/resource/hpo_genes.tsv.gz')
  
  sapply(1:nrow(sa2), function(i){
    hpos <- strsplit(sa2[i, HPO_TERMS], split = ',') %>% unlist
    genes <- paste(sort(hpo_dt[HPO_id %in% hpos, hgncSymbol]), collapse = ',')
    set(sa2, i, 'HPO_matching_genes', value = genes)
  }) %>% invisible()  # don't print result
  sa2 <- sa2[, .(RNA_ID, HPO_matching_genes)]
  
  fwrite(sa2, 
         file.path(snakemake@config$root, 'processed_data/sample_anno/genes_overlapping_HPO_terms.tsv'), 
         na = NA, sep = "\t", row.names = F, quote = F)
}

file.create(snakemake@output$done)

