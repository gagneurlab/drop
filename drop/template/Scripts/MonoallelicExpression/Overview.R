#'---
#' title: Monoallelic Expression
#' author:
#' wb:
#'  log:
#'    - snakemake: '`sm str(tmp_dir / "MAE" / "Overview.Rds")`'
#'  input:
#'    - allelic_counts: '`sm expand(cfg.getProcessedDataDir() + 
#'                          "/mae/allelic_counts/{mae_id}.csv.gz",
#'                          mae_id=cfg.MAE.getMaeAll())`'
#'    - results_obj: '`sm expand(cfg.getProcessedResultsDir() + 
#'                       "/mae/samples/{mae_id}_res.Rds", 
#'                       mae_id=cfg.MAE.getMaeAll())`'
#'    - results_tables: '`sm expand(cfg.getProcessedResultsDir() + 
#'                       "/mae/{dataset}/MAE_results_{annotation}.tsv", 
#'                       dataset=cfg.MAE.groups, annotation=cfg.getGeneVersions())`'
#'    - qc_matrix: '`sm expand(cfg.getProcessedResultsDir() + "/mae/{qc_group}/" +
#'                  "dna_rna_qc_matrix.Rds", qc_group=cfg.MAE.qcGroups)`'
#' output:
#'   html_document:
#'    code_folding: hide
#'    code_download: TRUE
#'---

#+ echo=F
saveRDS(snakemake, snakemake@log$snakemake)

suppressPackageStartupMessages({
  library(tMAE)
})


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

