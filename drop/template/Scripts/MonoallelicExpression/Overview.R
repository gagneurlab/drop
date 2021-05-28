#'---
#' title: Monoallelic Expression
#' author:
#' wb:
#'  log:
#'    - snakemake: '`sm str(tmp_dir / "MAE" / "Overview.Rds")`'
#'  params:
#'    - run: '`sm cfg.MAE.run`'
#'  input:
#'    - '`sm **aberrantSplicing_Overview_R_input(cfg)`'
#'  output:
#'   html_document:
#'    code_folding: hide
#'    code_download: TRUE
#'---

#+ include=FALSE
knitr::opts_chunk$set(eval = snakemake@params$run)

#+ eval=TRUE, echo=FALSE
saveRDS(snakemake, snakemake@log$snakemake)



#' ## Files
#' ### Allelic counts
#' Located in `r file.path(snakemake@config$root, 'processed_data/mae/allelic_counts/')`
#' 
#' ### Results tables of each sample
#' Located in `r file.path(snakemake@config$root, 'processed_results/mae/samples/')`
#' 
#' ### Aggregated results tables of each group
#' `r paste('* ', snakemake@input$results_tables, collapse = '\n')`  
#'
#' ### MAE Pipeline Output
#' [MAE Pipeline Output](`r "./Scripts_MAE_Datasets.html"`)
#' 

#' ## Analyze Individual Results
# Read the first results table
res_sample <- readRDS(snakemake@input$results_obj[[1]])

#+echo=F
library(tMAE)

if(is.null(res_sample$rare)){
  g1 <- plotMA4MAE(res_sample)
  g2 <- plotAllelicCounts(res_sample)
} else {
  g1 <- plotMA4MAE(res_sample, rare_column = 'rare')
  g2 <- plotAllelicCounts(res_sample, rare_column = 'rare')
}

#' ### MA plot: fold change vs RNA coverage
g1

#' ### Alternative vs Reference plot
g2

#' ## Quality Control: VCF-BAM Matching
#' 
#' [QC Overview](`r "./Scripts_QC_Datasets.html"`)
#' 
#' ### DNA-RNA matrix: 
#' `r paste('* ', snakemake@input$qc_matrix, collapse='\n')`  

