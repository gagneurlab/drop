#'---
#' title: Monoallelic Expression
#' author:
#' wb:
#'  log:
#'    - snakemake: '`sm str(tmp_dir / "MAE" / "Overview.Rds")`'
#'  params:
#'    - annotations: '`sm cfg.genome.getGeneVersions()`'
#'    - datasets: '`sm cfg.MAE.groups`'
#'    - qc_groups: '`sm cfg.MAE.qcGroups`'
#'    - htmlDir: '`sm config["htmlOutputPath"] + "/MonoallelicExpression"`'
#'    - resultsDir: '`sm cfg.getProcessedResultsDir() + "/mae"`'
#'  input:
#'    - functions: '`sm cfg.workDir / "Scripts/html_functions.R"`'
#'    - allelic_counts: '`sm expand(cfg.getProcessedDataDir() +
#'                          "/mae/allelic_counts/{mae_id}.csv.gz",
#'                          mae_id=cfg.MAE.getMaeAll())`'
#'    - results_obj: '`sm expand(cfg.getProcessedResultsDir() +
#'                       "/mae/samples/{mae_id}_res.Rds",
#'                       mae_id=cfg.MAE.getMaeAll())`'
#'    - results_tables: '`sm expand(cfg.getProcessedResultsDir() +
#'                       "/mae/{dataset}/MAE_results_{annotation}.tsv",
#'                       dataset=cfg.MAE.groups, annotation=cfg.genome.getGeneVersions())`'
#'    - qc_matrix: '`sm expand(cfg.getProcessedResultsDir() + "/mae/{qc_group}/" +
#'                  "dna_rna_qc_matrix.Rds", qc_group=cfg.MAE.qcGroups)`'
#' output:
#'   html_document:
#'    code_folding: show
#'    code_download: TRUE
#'---

#+ include=FALSE
saveRDS(snakemake, snakemake@log$snakemake)
source(snakemake@input$functions)

#+ eval=TRUE, echo=FALSE
# get parameters
datasets <- sort(snakemake@params$datasets)
annotations <- snakemake@params$annotations
htmlDir <- snakemake@params$htmlDir
resultsDir <- snakemake@params$resultsDir

results_links <- sapply(
  annotations, function(v) build_link_list(
    file_paths = file.path(htmlDir, paste0(datasets, '--', v, '_results.html')),
    captions = datasets
  )
)

table_links <- sapply(
  annotations, function(v) build_link_list(
    file_paths = file.path(resultsDir, paste0(datasets, '/MAE_results_', v, '.tsv')),
    captions = paste0(datasets)
  )
)

#'
#' **Datasets:** `r paste(datasets, collapse = ', ')`
#'
#' **Gene annotations:** `r paste(annotations, collapse = ', ')`
#'
#' ## MAE results
#' `r display_text(caption = 'Gene annotation version ', links = results_links)`
#'
#' ## Files
#' * [Allelic counts](`r file.path(snakemake@config$root, 'processed_data/mae/allelic_counts/')`)
#' * [Results data tables of each sample (.Rds)](`r file.path(snakemake@config$root, 'processed_results/mae/samples/')`)  

#'
#' `r display_text(caption = 'Significant MAE results tables ', links = table_links)`


#' ## Quality Control: VCF-BAM Matching
#+ eval=TRUE, echo=FALSE
qc_groups <- sort(snakemake@params$qc_groups)
qc_links <- build_link_list(
    file_paths = file.path(htmlDir, paste0('QC/', qc_groups, '.html')),
    captions = qc_groups
)

qc_matrix_links <- build_link_list(
    file_paths = file.path(snakemake@input$qc_matrix),
    captions = qc_groups
)

#' `r display_text(caption = 'QC Overview ', links = qc_links)`
#' `r display_text(caption = 'DNA-RNA matrix ', links = qc_matrix_links)`
#'

#+ eval=TRUE, echo=TRUE
#' ## Analyze Individual Results
# Read the first results table
res_sample <- readRDS(snakemake@input$results_obj[[1]])

#+echo=F
library(tMAE)

if(is.na(res_sample$rare)){
  g1 <- plotMA4MAE(res_sample)
  g2 <- plotAllelicCounts(res_sample)
} else {
  g1 <- plotMA4MAE(res_sample, rare_column = 'rare')
  g2 <- plotAllelicCounts(res_sample, rare_column = 'rare')
}

#' ### MA plot: fold change vs RNA coverage
#+echo=F
g1

#' ### Alternative vs Reference plot
#+echo=F
g2
