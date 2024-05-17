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
#'                       "/mae/allelic_counts/{mae_id}.csv.gz",
#'                        mae_id=cfg.MAE.getMaeAll())`'
#'    - results_sample: '`sm expand(cfg.getProcessedResultsDir() +
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
allelicRatioCutoff <- snakemake@config$mae$allelicRatioCutoff
padjCutoff <- snakemake@config$mae$padjCutoff

results_links <- sapply(annotations, function(x) build_link_list(
  file_paths = file.path(htmlDir, paste0(datasets, '--', x, '_results.html')),
  captions = datasets
)
)

table_links <- sapply(annotations, function(v) build_link_list(
  file_paths = file.path(resultsDir, paste0(datasets, '/MAE_results_', v, '.tsv')),
  captions = paste0(datasets)
)
)


#'
#' **Datasets:** `r paste(datasets, collapse = ', ')`
#'
#' **Gene annotations:** `r paste(annotations, collapse = ', ')`
#'
#' ## MAE results html
#' `r display_text(caption = 'Gene annotation ', links = results_links)`
#'
#' ## MAE files
#' * Allelic counts located under: `r file.path(snakemake@config$root, 'processed_data/mae/allelic_counts/')`
#' * Results for each sample located under: `r file.path(snakemake@config$root, 'processed_results/mae/samples/')`
#' * Aggregated results located under: `r file.path(snakemake@config$root, 'processed_results/mae/{DROP_group}/')`
#'

#' ## Quality Control: VCF-BAM matching html
#+ eval=TRUE, echo=FALSE
qc_groups <- sort(snakemake@params$qc_groups)
qc_links <- sapply(qc_groups, function(v) build_link_list(
  file_paths = file.path(htmlDir, paste0('QC/', v, '.html')),
  captions = paste0(qc_groups)
)
)

#' `r display_text(caption = 'QC Overview ', links = qc_links)`
#'

#' ## QC files
#' * QC matrix located under: `r file.path(snakemake@config$root, 'processed_results/mae/{DROP_group}/')`
#'

#' ## Analyze Individual Results
#+ echo=FALSE
# Read the first results table
res_sample <- readRDS(snakemake@input$results_sample[[1]])
sample <- unique(res_sample$ID)

library(tMAE)
library(ggplot2)
rare_column <- 'rare'
if(any(is.na(res_sample$rare))) rare_column <- NULL
#+ echo=TRUE

#' ### MA plot: fold change vs RNA coverage
plotMA4MAE(res_sample, rare_column = rare_column,
           padjCutoff = padjCutoff,
           allelicRatioCutoff = allelicRatioCutoff) + ggtitle(sample)

#' ### Alternative vs Reference plot
plotAllelicCounts(res_sample, rare_column = rare_column,
                  padjCutoff = padjCutoff,
                  allelicRatioCutoff = allelicRatioCutoff) + ggtitle(sample)
