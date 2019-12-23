#'---
#' title: Monoallelic Expression
#' author:
#' wb:
#'  params:
#'    - tmpdir: '`sm drop.getTmpDir()`'
#'    - mae_ids: '`sm parser.getMaeAll()`'
#'    - count_matrices: '`sm expand(parser.getProcDataDir() + 
#'                       "/mae/allelic_counts/{mae_id}.csv.gz",
#'                       mae_id=parser.getMaeAll())`'
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
#'    - qc_html: '`sm config["htmlOutputPath"] + "/Scripts_QC_Overview.html"`'
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
#' # All Samples
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
#' DNA-RNA matrix: 
#' 
#' `r paste('    *', snakemake@params$qc_matrix, collapse='\n')`  
#' 
#' [QC Overview](`r snakemake@params$html`)

