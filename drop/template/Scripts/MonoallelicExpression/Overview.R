#'---
#' title: Monoallelic Expression
#' author:
#' wb:
#'  log:
#'    - snakemake: '`sm str(tmp_dir / "MAE" / "Overview.Rds")`'
#'  params:
#'    - annotations: '`sm cfg.genome.getGeneVersions()`'
#'    - datasets: '`sm cfg.MAE.groups`'
#'    - htmlDir: '`sm config["htmlOutputPath"] + "/MAE"`'
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
#'    code_folding: hide
#'    code_download: TRUE
#'---

#+ include=FALSE

#+ eval=TRUE, echo=FALSE
saveRDS(snakemake, snakemake@log$snakemake)
source(snakemake@input$functions)

# get parameters
datasets <- sort(snakemake@params$datasets)
annotations <- snakemake@params$annotations
htmlDir <- snakemake@params$htmlDir

results_links <- sapply(
  annotations, function(v) build_link_list(
    file_paths = file.path(htmlDir, paste0(datasets, '--', v, '_results.html')),
    captions = datasets
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
#' * [Results tables of each sample](`r file.path(snakemake@config$root, 'processed_results/mae/samples/')`)
#' * [Aggregated results tables of each group](`r paste('* ', snakemake@input$results_tables, collapse = '\n')`)
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
