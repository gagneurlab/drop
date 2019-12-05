#'---
#' title: Analysis Example
#' author: salazar, mumichae, vyepez
#' wb:
#'  py:
#'    - |
#'     datasets = config["mae"]["groups"]
#'     annotations = list(config["geneAnnotation"].keys())
#'     mae_ids = parser.getMaeAll()
#'  params:
#'    - tmpdir: '`sm drop.getTmpDir()`'
#'    - mae_ids: '`sm mae_ids`'
#'    - count_matrices: '`sm expand(parser.getProcDataDir() + 
#'                       "/mae/allelic_counts/{mae_id}.csv.gz",
#'                       mae_id=mae_ids)`'
#'    - results_obj: '`sm expand(parser.getProcResultsDir() + 
#'                       "/mae/samples/{mae_id}_res.Rds", mae_id=mae_ids)`'
#'    - results_tables: '`sm expand(parser.getProcResultsDir() + 
#'                       "/mae/{dataset}/MAE_results_{annotation}.tsv", 
#'                       dataset=datasets, annotation=annotations)`'
#'    - html: '`sm config["htmlOutputPath"] + "/Scripts_MAE_Results_Overview.html"`'
#'    - qc_matrix: '`sm parser.getProcResultsDir() + "/mae/" +
#'                config["mae"]["qcGroup"] + "/dna_rna_qc_matrix.Rds"`'
#'    - qc_html: '`sm config["htmlOutputPath"] + "/Scripts_QC_DNA_RNA_matrix_plot.html"`'
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

all_files <- paste0("* ", snakemake@params$mae_ids, ": ",
       snakemake@params$count_matrices, collapse = "\n")
#'
#' # Monoallelic Expression
#' 
#'  `r all_files`
#' 
#' [MAE Pipeline Output](`r snakemake@params$html`)
#' 

#' # Analyze Individual Results
#' 
file <- snakemake@params$results_tables[[1]]
res <- fread(file)

file_location <- strsplit(file, "/")[[1]]
res_sample <- readRDS(snakemake@params$results_obj[[1]])

plotMA4MAE(res_sample)
plotAllelicCounts(res_sample)

#' # Quality Control: VCF-BAM Matching
#' 
#' DNA-RNA matrix: `r snakemake@params$matrix`  
#' 
#' [QC Overview](`r snakemake@params$html`)

