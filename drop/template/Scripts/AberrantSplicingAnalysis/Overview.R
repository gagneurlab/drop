#'---
#' title: Aberrant Splicing
#' author:
#' wb:
#'  py:
#'    - |
#'      datasets = config['aberrantSplicing']['groups']
#'  params:
#'    - tmpdir: '`sm drop.getTmpDir()`'
#'    - fds_files: '`sm expand(parser.getProcDataDir() +
#'                "/aberrant_splicing/datasets/savedObjects/{dataset}/" + 
#'                "fds-object.RDS", dataset=datasets)`'
#'    - result_tables: '`sm expand(parser.getProcDataDir() +
#'                    "/aberrant_splicing/results/{dataset}_results_per_junction.tsv",
#'                    dataset=datasets)`'
#'  input:
#'    - AS: '`sm drop.getTmpDir() + "/AS.done"`'
#' output:
#'   html_document:
#'    code_folding: hide
#'    code_download: TRUE
#'---

#+ echo=F
saveRDS(snakemake, file.path(snakemake@params$tmpdir, 
                             "AberrantSplicing_FRASER.snakemake"))
# snakemake <- readRDS(".drop/tmp/AberrantSplicing_FRASER.snakemake")

suppressPackageStartupMessages({
    library(FRASER)
    library(magrittr)
})

#' FRASER objects: `r paste(snakemake@params$fds_files)`  
#' 
#' FRASER results: `r paste(snakemake@params$result_tables)`  
#' 

#'
#' ## Analyze individual results
#' 

fds <- loadFraseRDataSet(file = snakemake@params$fds_files[[1]])
fds <- readRDS(snakemake@params$fds_files[[1]])
res <- fread(snakemake@params$result_tables[[1]])

sample <- samples(fds)[1]

#' Get a splice site and sample of interest
sample <- res[1, sampleID]
siteIndex <- 4

#' ## Volcano plot
#' setting basePlot = FALSE creates an interactive plot
#' that allows finding the junction(s) of interest
FRASER::plotVolcano(fds, sample, type = 'psi3', basePlot = TRUE)

#' ## Expression plot
FRASER::plotExpression(fds, type = 'psi3', site = siteIndex, basePlot = TRUE)

#' ## Expected vs observed PSI
FRASER::plotExpectedVsObservedPsi(fds, type = 'psi3', 
                                  idx = siteIndex, basePlot = TRUE)

