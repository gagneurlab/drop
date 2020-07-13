#'---
#' title: Aberrant Splicing
#' author:
#' wb:
#'  log:
#'    - snakemake: '`sm str(tmp_dir / "AS" / "Overview.Rds")`'
#'  input:
#'    - fds_files: '`sm expand(cfg.getProcessedDataDir() +
#'                "/aberrant_splicing/datasets/savedObjects/{dataset}/" + 
#'                "fds-object.RDS", dataset=cfg.AS.groups)`'
#'    - result_tables: '`sm expand(cfg.getProcessedDataDir() +
#'                    "/aberrant_splicing/results/{dataset}_results_per_junction.tsv",
#'                    dataset=cfg.AS.groups)`'
#' output:
#'   html_document:
#'    code_folding: show
#'    code_download: TRUE
#'---

#+ echo=F
saveRDS(snakemake, snakemake@log$snakemake)

suppressPackageStartupMessages({
  library(FRASER)
  library(magrittr)
})

# define functions
get_html_path <- function(datasets, htmlDir, fileName) {
  file_paths <- file.path(htmlDir, fileName)
  file_link <- paste0('\n* [', datasets ,'](', file_paths, 
                      '){target="_blank"}\n', collapse = ' ')
  file_link
}

display_text <- function(links) {
  paste0(links, collapse = '\n')
}

# get parameters
datasets <- sort(snakemake@config$aberrantSplicing$groups)

## start html

#'
#' **Datasets:** `r paste(datasets, collapse = ', ')`
#'
#' ## Summaries
#' ### Counts summary
#+ echo=FALSE
htmlDir <- './AberrantSplicing'
count_links <- get_html_path(datasets = datasets,
                             htmlDir = htmlDir, 
                             fileName = paste0(datasets, '_countSummary.html'))
#' 
#' `r display_text(count_links)`
#' 
#' ### FRASER summary
#+ echo=FALSE
fraser_links <- get_html_path(datasets = datasets,
                              htmlDir = htmlDir, 
                              fileName = paste0(datasets, '_summary.html'))
#' 
#' `r display_text(fraser_links)`

#' ## Files
#' ### FRASER datasets (fds)
#' `r paste('* ', snakemake@input$fds_files, collapse = '\n')`  
#' 
#' ### Results tables
#' `r paste('* ', snakemake@input$result_tables, collapse = '\n')`  

#'
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
#' setting basePlot = FALSE creates an interactive plot
#' that allows finding the junction(s) of interest
FRASER::plotVolcano(fds, sample, type = 'psi3', basePlot = TRUE)

#' ### Expression plot
FRASER::plotExpression(fds, type = 'psi3', site = siteIndex, basePlot = TRUE)

#' ### Expected vs observed PSI
FRASER::plotExpectedVsObservedPsi(fds, type = 'psi3', 
                                  idx = siteIndex, basePlot = TRUE)

