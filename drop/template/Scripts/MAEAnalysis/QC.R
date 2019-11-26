#'---
#' title: OUTRIDER results
#' author: mumichae
#' wb:
#'  input:
#'    - matrix: '`sm parser.getProcResultsDir() + "/mae/" + config["mae"]["qcGroup"] + "/dna_rna_qc_matrix.Rds"`'
#'    - html: '`sm config["htmlOutputPath"] + "/Scripts_QC_DNA_RNA_matrix_plot.html"`'
#' output:
#'   html_document:
#'    code_folding: hide
#'    code_download: TRUE
#'---

print(getwd())
#+ echo=F
saveRDS(snakemake, '.tmp/mae_qc.snakemake')
# snakemake <- readRDS('.tmp/mae_qc.snakemake')

#' `r snakemake@input$matrix`  
#' `r snakemake@input$html`
