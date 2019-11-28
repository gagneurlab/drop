#'---
#' title: Analysis Example
#' author: salazar
#' wb:
#'  py:
#'    - |
#'     datasets = config["mae"]["groups"]
#'     annotations = list(config["geneAnnotation"].keys())
#'     mae_ids = parser.getMaeAll()
#'  params:
#'    - tmpdir: '`sm drop.getTmpDir()`'
#'  input:
#'    - count_matrices: '`sm expand(parser.getProcDataDir() + 
#'                       "/mae/allelic_counts/{mae_id}.csv.gz",
#'                       mae_id=mae_ids)`'
#'    - results_tables: '`sm expand(parser.getProcResultsDir() + 
#'                       "/mae/{dataset}/MAE_results_{annotation}.tsv", 
#'                       dataset=datasets, annotation=annotations)`'
#'    - html: '`sm config["htmlOutputPath"] + "/Scripts_MAE_Results_Overview.html"`'
#' output:
#'   html_document:
#'    code_folding: hide
#'    code_download: TRUE
#'---

#+ echo=F
saveRDS(snakemake, file.path(snakemake@params$tmpdir, "MAE_analysis.snakemake"))
# snakemake <- readRDS(".drop/tmp/MAE_analysis.snakemake")

#' `r snakemake@input$matrix`  
#' `r snakemake@input$html`


library(tMAE)
#file <- results_tables[[1]]
#res <- fread(file)
#sample <- res[1, MAE_ID]

#' Load the file of interest
#file_location <- strsplit(file, "/")[[1]]
#res_sample <- readRDS(paste0(paste(file_location[1:eval(length(file_location)-2)], collapse = "/"), "/samples/", sample, "_res.Rds"))

# Plots
#plotMA4MAE(res_sample, rare_column = 'rare')
#plotAllelicCounts(res_sample, rare_column = 'rare')
