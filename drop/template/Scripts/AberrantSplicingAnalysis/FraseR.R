#'---
#' title: Analysis Example
#' author: salazar, mumichae
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


#'    - fds_files: '`sm expand(parser.getProcDataDir() + "/aberrant_splicing/datasets/savedObjects/{dataset}/pajdBinomial_psiSite.h5", dataset=datasets)`'
#'    - result_tables: '`sm expand(parser.getProcDataDir() + "/aberrant_splicing/results/{dataset}_results.tsv", dataset=datasets)`'

#' FraseR objects: `r paste(snakemake@input$fds_files)`  
#' FraseR results: `r paste(snakemake@input$result_tables)`  


