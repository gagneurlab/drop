#'---
#' title: Sample Annotation Overview
#' author:
#' wb:
#'  log:
#'   - snakemake: '`sm str(tmp_dir / "SampleAnnotation.Rds")`'
#'  params:
#'   - export_dir: '`sm cfg.getProcessedResultsDir() + "/exported_counts"`'
#'   - groups: '`sm cfg.exportCounts.getExportGroups()`'
#'   - hpoFile: '`sm cfg.get("hpoFile")`'
#'  input: 
#'   - sampleAnnotation: '`sm config["sampleAnnotation"]`'
#'  output:
#'   - export: '`sm touch(cfg.getProcessedResultsDir() + "/exported_counts/sample_anno.done")`'
#'   - done: '`sm touch(cfg.getProcessedDataDir() + "/sample_anno/sample_anno.done")`'
#' output:
#'   html_document:
#'    code_folding: hide
#'    code_download: TRUE
#'---

#+echo=F
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
if('DNA_VCF_FILE' %in% colnames(sa)){
  sa[, aux1 := file.exists(DNA_VCF_FILE) | is.na(DNA_ID)]
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
#+echo=F
if(!is.null(sa$HPO_TERMS)){
  if(!all(is.na(sa$HPO_TERMS)) & ! all(sa$HPO_TERMS == '')){
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
  
  fwrite(sa2, 
         file.path(snakemake@config$root, 
                   'processed_data/sample_anno/genes_overlapping_HPO_terms.tsv'), 
         na = NA, sep = "\t", row.names = F, quote = F)
  }
}
sa[, DROP_GROUP := gsub(' ', '', DROP_GROUP)]
if(!is.null(sa$ICD_10))
  sa[, ICD_10 := toupper(ICD_10)]

list_groups <- strsplit(sa$DROP_GROUP, split = ',')
drop_groups <- snakemake@params$groups # list_groups %>% unlist %>% unique

genomeAssembly <- snakemake@config$genomeAssembly
geneAnnotation <- snakemake@config$exportCounts$geneAnnotations

sa[STRAND == 'no', STRAND_SPECIFIC := FALSE]
sa[STRAND %in% c('yes', 'reverse'), STRAND_SPECIFIC := TRUE]
sa[, STRAND := NULL]

# export processed sample annotation and DESCRIPTION file
for(dataset in drop_groups){
  cols <- intersect(c('RNA_ID', 'INDIVIDUAL_ID', 'TISSUE', 'SEX', 'AFFECTED', 
                      'ICD_10', 'PAIRED_END', 'STRAND_SPECIFIC'),
                    colnames(sa))
  sa_sub <- sa[sapply(list_groups, function(x) dataset %in% x), cols, with = F]
  path <- file.path(snakemake@params$export_dir, paste(dataset, genomeAssembly, geneAnnotation, sep = '--'))
  dir.create(path)
  fwrite(sa_sub, file = file.path(path, 'sampleAnnotation.tsv'), 
         quote = FALSE, row.names = FALSE, sep = '\t')
  
  # DESCRIPTION file
  v0 <- 'Title: # Add a title'
  v1 <- paste0('Number of samples: ', nrow(sa_sub))
  v2 <- ifelse(!is.null(sa_sub$TISSUE), 
               ifelse(uniqueN(sa_sub$TISSUE) > 1, 
                      'More than 1 tissue in dataset is not allowed!',
                      paste0('Tissue: ', unique(sa_sub$TISSUE))),
               'No tissue information')
  v3 <- 'Organism: Homo sapiens'
  v4 <- paste0('Genome assembly: ', genomeAssembly)
  v5 <- paste0('Gene annotation: ', geneAnnotation)
  v6 <- ifelse(!is.null(sa_sub$ICD_10), 
               paste0('Disease (ICD-10: N): ', 
                      paste(unite(data.table(table(sa_sub$ICD_10)), col = 'aux', 'V1', 'N', sep = ': ')$aux, collapse = ', ' )),
               'No disease information')
  v7 <- ifelse(uniqueN(sa_sub$STRAND_SPECIFIC) > 1, 
               'All samples should be either strand- or non-strand-specific!',
               paste0('Strand specific: ', unique(sa_sub$STRAND_SPECIFIC)))
  v8 <- ifelse(uniqueN(sa_sub$PAIRED_END) > 1, 
               'All samples should be either single end or paired end!',
               paste0('Paired end: ', unique(sa_sub$PAIRED_END)))
  v9 <- 'Cite as: RNA-Seq count tables were taken from: # add your citation(s)'
  v10 <- 'Dataset contact: # Use format Name Last_Name, <email address>'
  v11 <- 'Comments: # add any comments, if needed, otherwise remove'
  
  writeLines(c(v0, v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11), file.path(path, 'DESCRIPTION.txt'))
}

