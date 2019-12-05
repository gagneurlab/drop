#'---
#' title: Analysis Example
#' author: salazar, mumichae, vyepez
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
#'                    "/aberrant_splicing/results/{dataset}_results.tsv",
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
    library(FraseR)
    library(magrittr)
})

#' FraseR objects: `r paste(snakemake@params$fds_files)`  
#' 
#' FraseR results: `r paste(snakemake@params$result_tables)`  
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
#' Hover over the plot and find your splice site(s) of interest
FraseR::plotVolcano(fds, sample, type = 'psi3')

#' ## Gene expression plot
#' Hover over the plot and find your sample(s) of interest
FraseR::plotExpression(fds, type = 'psi3', site = siteIndex)

FraseR::plotExpectedVsObservedPsi(fds, type = 'psi3', idx = siteIndex)

