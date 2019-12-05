#'---
#' title: OUTRIDER results
#' author: mumichae
#' wb:
#'  params:
#'    - tmpdir: '`sm drop.getTmpDir()`'
#'    - matrix: '`sm parser.getProcResultsDir() + "/mae/" +
#'                config["mae"]["qcGroup"] + "/dna_rna_qc_matrix.Rds"`'
#'    - html: '`sm config["htmlOutputPath"] + "/Scripts_QC_DNA_RNA_matrix_plot.html"`'
#'  input:
#'    - MAE: '`sm drop.getTmpDir() + "/MAE.done"`'
#' output:
#'   html_document:
#'    code_folding: hide
#'    code_download: TRUE
#'---

print(getwd())
#+ echo=F
saveRDS(snakemake, file.path(snakemake@params$tmpdir, 'mae_qc.snakemake'))
# snakemake <- readRDS('.tmp/mae_qc.snakemake')

#' `r snakemake@params$matrix`  
#' `r snakemake@params$html`
