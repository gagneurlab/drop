#'---
#' title: Aberrant Splicing
#' author:
#' wb:
#'  log:
#'    - snakemake: '`sm str(tmp_dir / "AS" / "Overview.Rds")`'
#'  params:
#'    - annotations: '`sm cfg.genome.getGeneVersions()`'
#'    - datasets: '`sm cfg.AS.groups`'
#'    - htmlDir: '`sm config["htmlOutputPath"] + "/AberrantSplicing"`'
#'  input:
#'    - functions: '`sm cfg.workDir / "Scripts/html_functions.R"`'
#'    - fds_files: '`sm expand(cfg.getProcessedResultsDir() +
#'                "/aberrant_splicing/datasets/savedObjects/{dataset}--{annotation}/" +
#'                "fds-object.RDS", dataset=cfg.AS.groups, annotation=cfg.genome.getGeneVersions())`'
#'    - result_tables: '`sm expand(cfg.getProcessedResultsDir() +
#'                    "/aberrant_splicing/results/{annotation}/fraser/{dataset}/results_per_junction.tsv",
#'                    dataset=cfg.AS.groups, annotation=cfg.genome.getGeneVersions())`'
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

count_links <- build_link_list(
  file_paths = file.path(htmlDir, paste0(datasets, '_countSummary.html')),
  captions = datasets
)

results_links <- sapply(
  annotations, function(v) build_link_list(
    file_paths = file.path(htmlDir, paste0(datasets, '--', v, '_summary.html')),
    captions = datasets
  )
)

fds_links <- build_link_list(snakemake@input$fds_files)
results_tables <- build_link_list(snakemake@input$result_tables)

## start html

#'
#' **Datasets:** `r paste(datasets, collapse = ', ')`
#'
#' **Gene annotations:** `r paste(annotations, collapse = ', ')`
#'
#' ## Summaries
#' ### Counts summary
#' `r display_text(links = count_links)`
#'
#' ### FRASER summary
#' `r display_text(caption = 'Gene annotation version ', links = results_links)`
#'
#' ## Files
#' `r display_text(caption = 'FRASER datasets (fds)', links = fds_links)`
#' `r display_text(caption = 'Results tables', links = results_tables)`
#'

#+ echo=FALSE
library(FRASER)
library(magrittr)

#' ## Analyze individual results
# Read the first fds object and results table
fds <- loadFraserDataSet(file = snakemake@input$fds_files[[1]])
res <- fread(snakemake@input$result_tables[[1]])

#' Display the results table of the first dataset
#+ echo=FALSE
DT::datatable(res, filter = 'top')

#' Get a splice site and sample of interest. Outliers are in red.
#+ echo=TRUE
sample <- res[1, sampleID]
siteIndex <- 4

#' ### Volcano plot
# set basePlot to FALSE to create an interactive plot
FRASER::plotVolcano(fds, sample, type = 'psi3', basePlot = TRUE)

#' ### Expression plot
FRASER::plotExpression(fds, type = 'psi3', site = siteIndex, basePlot = TRUE)

#' ### Expected vs observed PSI (or theta)
FRASER::plotExpectedVsObservedPsi(fds, type = 'psi3',
                                  idx = siteIndex, basePlot = TRUE)
