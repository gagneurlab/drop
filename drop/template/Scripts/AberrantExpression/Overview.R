#'---
#' title: Aberrant Expression
#' author:
#' wb:
#'  log:
#'    - snakemake: '`sm str(tmp_dir / "AE" / "Overview.Rds")`'
#'  params:
#'    - annotations: '`sm cfg.genome.getGeneVersions()`'
#'    - datasets: '`sm cfg.AE.groups`'
#'    - htmlDir: '`sm config["htmlOutputPath"] + "/AberrantExpression"`'
#'  input:
#'    - functions: '`sm cfg.workDir / "Scripts/html_functions.R"`'
#'    - odsFiles: '`sm expand(cfg.getProcessedResultsDir() +
#'                  "/aberrant_expression/{annotation}/outrider/{dataset}/ods.Rds",
#'                  annotation=cfg.genome.getGeneVersions(), dataset=cfg.AE.groups)`'
#'    - resultTables: '`sm expand(cfg.getProcessedResultsDir() +
#'                      "/aberrant_expression/{annotation}/outrider/" +
#'                      "{dataset}/OUTRIDER_results.tsv",
#'                      annotation=cfg.genome.getGeneVersions(), dataset=cfg.AE.groups)`'
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

count_links <- sapply(
  annotations, function(x) build_link_list(
    file_paths = file.path(htmlDir, "Counting", x, paste0('Summary_', datasets, '.html')),
    captions = datasets)
)

results_links <- sapply(
 annotations, function(x) build_link_list(
    file_paths = file.path(htmlDir, "Outrider", x, paste0('Summary_', datasets, '.html')),
    captions = datasets)
)

ods_links <- build_link_list(snakemake@input$odsFiles)
results_tables <- build_link_list(snakemake@input$resultTables)

## start html

#'
#' **Datasets:** `r paste(datasets, collapse = ', ')`
#'
#' **Gene annotations:** `r paste(annotations, collapse = ', ')`
#'
#' ## Summaries
#'
#' ### Counts summary
#'
#' `r display_text(caption = 'Gene annotation version ', links = count_links)`
#'
#' ### OUTRIDER summary
#'
#' `r display_text(caption = 'Gene annotation version ', links = results_links)`
#'
#' ## Files
#' `r display_text(caption = 'OUTRIDER datasets (ods)', links = ods_links)`
#' `r display_text(caption = 'Results tables', links = results_tables)`
#'

#' ## Analyze Individual Results
#+ echo=FALSE
library(OUTRIDER)

# Read the first ods object and results table
ods <- readRDS(snakemake@input$odsFiles[[1]])
res <- fread(snakemake@input$resultTables[[1]])

#' Display the results table of the first dataset
#+ echo=FALSE
if(nrow(res) > 0){
DT::datatable(head(res, 100), caption = 'OUTRIDER results (up to 100 rows shown)',
                options=list(scrollX=TRUE), filter = 'top')
  
} else print("no significant results")

#' Choose a random gene and sample to plot. Outliers are in red.
#+ echo=TRUE
gene <- ifelse(nrow(res)>1, res[1, geneID], rownames(ods)[1])
sample <- ifelse(nrow(res)>1, res[1, sampleID], colnames(ods)[1])

#' ### Volcano plot
#' setting basePlot = FALSE creates an interactive plot
#' that allows finding the gene(s) of interest
OUTRIDER::plotVolcano(ods, sample, basePlot = TRUE,
                      zScoreCutoff = snakemake@config$aberrantExpression$zScoreCutoff,
                      padjCutoff = snakemake@config$aberrantExpression$padjCutoff)

#' ### Gene expression plot (normalized counts)
OUTRIDER::plotExpressionRank(ods, gene, basePlot = TRUE,
                      zScoreCutoff = snakemake@config$aberrantExpression$zScoreCutoff,
                      padjCutoff = snakemake@config$aberrantExpression$padjCutoff)

#' ### Expected vs observed counts
OUTRIDER::plotExpectedVsObservedCounts(ods, gene, basePlot = TRUE)
