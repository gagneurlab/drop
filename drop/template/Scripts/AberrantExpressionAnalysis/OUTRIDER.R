#'---
#' title: OUTRIDER results
#' author: mumichae, vyepez
#' wb:
#'  py:
#'   - |
#'    annotations = list(config["geneAnnotation"].keys())
#'    datasets = config["aberrantExpression"]["groups"]
#'  params:
#'    - tmpdir: '`sm drop.getTmpDir()`'
#'  input:
#'    - ods_files: '`sm expand(parser.getProcResultsDir() + "/aberrant_expression/{annotation}/outrider/{dataset}/ods.Rds", annotation=annotations, dataset=datasets)`'
#'    - result_tables: '`sm expand(parser.getProcResultsDir() + "/aberrant_expression/{annotation}/outrider/{dataset}/OUTRIDER_results.tsv", annotation=annotations, dataset=datasets)`'
#'    - html: '`sm expand(config["htmlOutputPath"] + "/AberrantExpression/Outrider/{annotation}/Summary_{dataset}.html", annotation=annotations, dataset=datasets)`'
#' output:
#'   html_document:
#'    code_folding: hide
#'    code_download: TRUE
#'---

#+ echo=F
saveRDS(snakemake, file.path(snakemake@params$tmpdir, 
                             "AberrantExpression_OUTRIDER.snakemake"))
# snakemake <- readRDS(".drop/tmp/AberrantExpression_OUTRIDER.snakemake")

groups <- names(snakemake@config$outrider_all)
anno_version <- names(snakemake@config$geneAnnotation)
html_file_dir <- gsub(snakemake@config$htmlOutputPath, '.', dirname(dirname(snakemake@input$html)))
summaries_titles <- sapply(anno_version, function(v) {
  paste0('[', groups ,'](', html_file_dir, '/', v, '/Summary_', groups, '.html){target="_blank"}', collapse = ' ')
}
)

#' ## OUTRIDER Results
#' The summaries can be found here:  
#' `r paste0('Gene annotation version ', names(summaries_titles), ': ', summaries_titles, collapse = '\n')`
#' 
#' Links to the OUTRIDER output and results files:
#' `r paste(snakemake@input$ods_files, collapse = '\n')`  
#' `r paste(snakemake@input$result_tables, collapse = '\n')`
#' 

#' ## Analyze individual results
#' ### Read outrider object and results
library(OUTRIDER)
#ods <- readRDS(snakemake@input$ods_files[[1]])
#res <- fread(snakemake@input$results_tables[[1]])
 
#' Get a gene and sample of interest
#gene <- res[1, geneID]
#sample <- res[1, sampleID]

#' Example of a volcano plot
#plotVolcano(ods, sample)  # scroll over the plot and find your gene(s) of interest

#' Gene expression plot
#plotExpressionRank(ods, gene)  # scroll over the plot and find your sample(s) of interest


