#'---
#' title: Aberrant Splicing
#' author:
#' wb:
#'  log:
#'    - snakemake: '`sm str(tmp_dir / "AS" / "Overview.Rds")`'
#'  params:
#'    - run: '`sm cfg.AS.run`'
#'    - annotations: '`sm cfg.genome.getGeneVersions()`'
#'    - datasets: '`sm cfg.AS.groups`'
#'    - htmlDir: '`sm config["htmlOutputPath"] + "/AberrantSplicing"`'
#'  input:
#'    - '`sm **aberrantSplicing_Overview_R_input(cfg)`'
#' output:
#'   html_document:
#'    code_folding: show
#'    code_download: TRUE
#'---


#+ include=FALSE
knitr::opts_chunk$set(eval = snakemake@params$run)

#+ echo=TRUE, eval = TRUE
saveRDS(snakemake, snakemake@log$snakemake)
print(snakemake@input)


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
datasets <- sort(snakemake@params$datasets)
annotations <- snakemake@params$annotations
htmlDir <- snakemake@params$htmlDir

count_links <- get_html_path(datasets = datasets,
                             htmlDir = htmlDir, 
                             fileName = paste0(datasets, '_countSummary.html'))

datasets_annotations <- as.character(outer(datasets, annotations, FUN = paste, sep = '--'))
fraser_links <- get_html_path(datasets = datasets_annotations,
                              htmlDir = htmlDir, 
                              fileName = paste0(datasets_annotations, '_summary.html'))
## start html

#'
#' **Datasets:** `r paste(datasets, collapse = ', ')`
#'
#' **Gene annotations:** `r paste(annotations, collapse = ', ')`
#'
#' ## Summaries
#' ### Counts summary
#+ echo=FALSE
# htmlDir <- './AberrantSplicing'
#' 
#' `r display_text(count_links)`
#' 
#' ### FRASER summary
#' 
#' `r display_text(fraser_links)`
#' 
#' ## Files
#' ### FRASER datasets (fds)
#' `r paste('* ', snakemake@input$fds_files, collapse = '\n')`  
#' 
#' ### Results tables
#' `r paste('* ', snakemake@input$result_tables, collapse = '\n')`  
#'
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
