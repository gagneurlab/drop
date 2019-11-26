#'---
#' title: Analysis Example
#' author: salazar, mumichae, vyepez
#' wb:
#'  py:
#'   - |
#'     datasets = config['aberrantSplicing']['groups']
#'  input:
#' output:
#'   html_document:
#'    code_folding: hide
#'    code_download: TRUE
#'---

#+ echo=F
saveRDS(snakemake, ".tmp/AberrantSplicing_Analysis.snakemake")
# snakemake <- readRDS(".tmp/AberrantSplicing_Analysis.snakemake")


#'  - fds_files: '`sm expand(parser.getProcDataDir() + "/aberrant_splicing/datasets/savedObjects/{dataset}/pajdBinomial_psiSite.h5", dataset=datasets)`'
#'  - result_tables: '`sm expand(parser.getProcDataDir() + "/aberrant_splicing/results/{dataset}_results.tsv", dataset=datasets)`'

#' FraseR objects: `r paste(snakemake@input$fds_files)`  
#' FraseR results: `r paste(snakemake@input$result_tables)`  


#' ## Analyze individual results
#' ### Read outrider object and results
library(FraseR)
library(magrittr)
file <- '/s/project/fraser/snakemake_pipeline/Data/paperPipeline/datasets/savedObjects/Adrenal_Gland__PCA/pajdBetaBinomial_psi3.h5'
fds <- loadFraseRDataSet(file = file)
# fds <- readRDS(fds_files[[1]])
res <- fread(results_tables[[1]])

#' check all samples
samples(fds) %>% sort
sample <- samples(fds)[1]

#' Get a gene and sample of interest
gene <- res[1, geneID]
sample <- res[1, sampleID]

#' Example of a volcano plot
plotVolcano(fds, sample, type = 'psi3')  # scroll over the plot and find your gene(s) of interest

#' Gene expression plot
plotExpression(fds, type = 'psi3', site = 42850)  # scroll over the plot and find your sample(s) of interest

plotExpectedVsObservedPsi(fds, type = 'psi3', idx = 42850)
