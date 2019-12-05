#'---
#' title: Overview
#' author: mumichae, vyepez
#' wb:
#'  py:
#'   - |
#'    annotations = list(config["geneAnnotation"].keys())
#'    datasets = config["aberrantExpression"]["groups"]
#'  params:
#'    - annotations: '`sm annotations`'
#'    - datasets: '`sm datasets`'
#'    - tmpdir: '`sm drop.getTmpDir()`'
#'    - htmlDir: '`sm config["htmlOutputPath"] + "/AberrantExpression"`'
#'    - odsFiles: '`sm expand(parser.getProcResultsDir() +
#'                  "/aberrant_expression/{annotation}/outrider/{dataset}/ods.Rds",
#'                  annotation=annotations, dataset=datasets)`'
#'    - resultTables: '`sm expand(parser.getProcResultsDir() +
#'                      "/aberrant_expression/{annotation}/outrider/" +
#'                      "{dataset}/OUTRIDER_results.tsv",
#'                      annotation=annotations, dataset=datasets)`'
#'  input:
#'    - AE: '`sm drop.getTmpDir() + "/AE.done"`'
#' output:
#'   html_document:
#'    code_folding: hide
#'    code_download: TRUE
#'---

#+ echo=F
saveRDS(snakemake, file.path(snakemake@params$tmpdir, 
                             "AberrantExpression_OUTRIDER.snakemake"))
# snakemake <- readRDS(".drop/tmp/AberrantExpression_OUTRIDER.snakemake")

suppressPackageStartupMessages({
    library(OUTRIDER)
})

# define functions
get_html_path <- function(annotationVersion, datasets, htmlDir, fileName) {
    file_paths <- file.path(htmlDir, annotationVersion, fileName)
    file_link <- paste0('\n* [', datasets ,'](', file_paths, 
                        '){target="_blank"}\n', collapse = ' ')
    file_link
}

display_text <- function(caption, links) {
    paste0('**', caption, '**', names(links), '\n', links, collapse = '\n')
}

# get parameters
datasets <- snakemake@params$datasets
annotations <- snakemake@params$annotations
htmlDir <- snakemake@params$htmlDir

## start html

#'
#' **Datasets:** `r paste(datasets, collapse = ', ')`
#' 
#' **Gene annotations:** `r paste(annotations, collapse = ', ')`
#'
#' # Count Summaries
#+ echo=FALSE
count_links <- sapply(annotations, get_html_path, 
                      datasets = datasets,
                      htmlDir = file.path(htmlDir, "Counting"), 
                      fileName = paste0('Summary_', datasets, '.html'))
#' 
#' `r display_text(caption = 'Gene annotation version ', count_links)`
#' 
#' # OUTRIDER Results
#+ echo=FALSE
outrider_links <- sapply(annotations, get_html_path, 
                      datasets = datasets,
                      htmlDir = file.path(htmlDir, "Outrider"), 
                      fileName = paste0('Summary_', datasets, '.html'))
#' 
#' 
#' `r display_text(caption = 'Gene annotation version ', outrider_links)`
#' 
#' 
#' OUTRIDER dataset (ods) files 
#' 
#' `r paste('* ', snakemake@params$odsFiles, collapse = '\n')`
#' 
#' 
#' Results tables
#' 
#'  `r paste('* ', snakemake@params$resultTables, collapse = '\n')`
#' 
#'   
#' # Analyze Individual Results
ods <- readRDS(snakemake@params$odsFiles[[1]])
res <- fread(snakemake@params$resultTables[[1]])

DT::datatable(res)

#' Choose a random gene and sample
#+ echo=TRUE
gene <- res[1, geneID]
sample <- res[1, sampleID]

#' ## Volcano plot
#' Hover over the plot and find your gene(s) of interest
OUTRIDER::plotVolcano(ods, sample)

#' ## Gene expression plot
#' Hover over the plot and find your sample(s) of interest
OUTRIDER::plotExpressionRank(ods, gene)

#' ## Expected vs observed counts
#' 
OUTRIDER::plotExpectedVsObservedCounts(ods, gene)

