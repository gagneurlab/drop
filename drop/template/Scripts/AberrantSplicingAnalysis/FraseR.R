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
#'                "/aberrant_splicing/datasets/savedObjects/{dataset}/pajdBetaBinomial_psiSite.h5",
#'                dataset=datasets)`'
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

#' FraseR objects: `r paste(snakemake@params$fds_files)`  
#' FraseR results: `r paste(snakemake@params$result_tables)`  


#' ## Analyze individual results
#' ### Read outrider object and results
library(FraseR)
library(magrittr)
# fds <- loadFraseRDataSet(file = snakemake@params$fds_files[1])
# fds <- readRDS(fds_files[[1]])
#res <- fread(results_tables[[1]])

#' check all samples
#samples(fds) %>% sort
#sample <- samples(fds)[1]

#' Get a gene and sample of interest
#gene <- res[1, geneID]
#sample <- res[1, sampleID]

#' Example of a volcano plot
#plotVolcano(fds, sample, type = 'psi3')  # scroll over the plot and find your gene(s) of interest

#' Gene expression plot
#plotExpression(fds, type = 'psi3', site = 42850)  # scroll over the plot and find your sample(s) of interest

#plotExpectedVsObservedPsi(fds, type = 'psi3', idx = 42850)
