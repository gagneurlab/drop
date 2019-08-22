#'---
#' title: Genetic Diagnosis Overview
#' author: mumichae, salazar
#' wb:
#'  input:
#'  - abExp: '`sm config["htmlOutputPath"] + "/Scripts_Outrider_Overview.html"  `'
#'  - mae:  '`sm config["htmlOutputPath"] + "/Scripts_MAE_Results_MAE.html"  `' 
#' output:
#'   html_document:
#'    code_folding: hide
#'    code_download: TRUE
#'---

### Index script
#' INDEX 
1
#gene_annotation_names <- names(snakemake@config$GENE_ANNOTATION)
#summarie_titles <- paste(gene_annotation_names, names(snakemake@config$outrider_filtered))
#summaries <- paste('[', summarie_titles ,'](', gsub(snakemake@config$htmlOutputPath, ".", snakemake@input$summaries), ')', sep = '')
#summaries <- paste(summaries, sep = '\n')

# 
# #'  - abb_expr: '`sm config["htmlOutputPath"] + "/aberrant-expression_index.html" `'
# #'  - abb_splicing: '`sm config["htmlOutputPath"] + "/aberrant-splicing_index.html" `'
# #'  - mae: '`sm config["htmlOutputPath"] + "/mae_index.html" `' #' output:

# #'  - abSpl: '`sm expand(parser.getProcDataDir() + "/aberrant_splicing/FraseR/{dataset}_results.html", dataset=config["DATASETS"] )`'


#### make it more automated!!!!! ---> automatically get all indexes that we have! 