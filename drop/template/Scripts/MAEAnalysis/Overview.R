#'---
#' title: Monoallelic Expression
#' author:
#' wb:
#'  params:
#'    - tmpdir: '`sm drop.getTmpDir()`'
#'    - mae_ids: '`sm parser.getMaeAll()`'
#'    - allelic_counts: '`sm expand(parser.getProcDataDir() + 
#'                          "/mae/allelic_counts/{mae_id}.csv.gz",
#'                          mae_id=parser.getMaeAll())`'
#'    - results_obj: '`sm expand(parser.getProcResultsDir() + 
#'                       "/mae/samples/{mae_id}_res.Rds", 
#'                       mae_id=parser.getMaeAll())`'
#'    - results_tables: '`sm expand(parser.getProcResultsDir() + 
#'                       "/mae/{dataset}/MAE_results_{annotation}.tsv", 
#'                       dataset=config["mae"]["groups"],
#'                       annotation=list(config["geneAnnotation"].keys()))`'
#'    - html: '`sm config["htmlOutputPath"] + "/Scripts_MAE_Results_Overview.html"`'
#'    - qc_matrix: '`sm expand(parser.getProcResultsDir() + "/mae/{qc_group}/" +
#'                  "dna_rna_qc_matrix.Rds", qc_group=config["mae"]["qcGroups"])`'
#'  input:
#'    - MAE: '`sm drop.getTmpDir() + "/MAE.done"`'
#' output:
#'   html_document:
#'    code_folding: hide
#'    code_download: TRUE
#'---

#+ echo=F
saveRDS(snakemake, file.path(snakemake@params$tmpdir, "MAE_analysis.snakemake"))
# snakemake <- readRDS(".drop/tmp/MAE_analysis.snakemake")

suppressPackageStartupMessages({
  library(tMAE)
})


#' ## Files
#' ### Allelic counts
#' Located in `r file.path(snakemake@config$root, 'processed_data/mae/allelic_counts/')`
#' 
#' ### Results tables of each sample
#' Located in `r file.path(snakemake@config$root, 'processed_results/mae/samples/')`
#' 
#' ### Aggregated results tables of each group
#' `r paste('* ', snakemake@params$results_tables, collapse = '\n')`  
#'
#' ### MAE Pipeline Output
#' [MAE Pipeline Output](`r "./Scripts_MAE_Datasets.html"`)
#' 

#' ## Analyze Individual Results
# Read the first results table
res_sample <- readRDS(snakemake@params$results_obj[[1]])

#+echo=F
if(is.null(res_sample$rare)){
  g1 <- plotMA4MAE(res_sample)
  g2 <- plotAllelicCounts(res_sample)
} else {
  g1 <- plotMA4MAE(res_sample, rare_column = 'rare')
  g2 <- plotAllelicCounts(res_sample, rare_column = 'rare')
}

#' ### MA plot: fold change vs RNA coverage
g1

#' ### Alternative vs Reference plot
g2

#' ## Quality Control: VCF-BAM Matching
#' 
#' [QC Overview](`r "./Scripts_QC_Datasets.html"`)
#' 
#' ### DNA-RNA matrix: 
#' `r paste('* ', snakemake@params$qc_matrix, collapse='\n')`  

