#'---
#' title: Results Overview
#' author: mumichae
#' wb:
#'  input:
#'  - variants: '`sm variants( "../variant-annotation-pipeline/Output/html/index.html")`'
#'  - abb_expr: '`sm aberrantExp( "../aberrant-expression-pipeline/Output/html/index.html")`'
#'  - abb_splicing: '`sm aberrantSplicing( "../aberrant-splicing-pipeline/Output/html/index.html")`'
#'  - mae: '`sm mae( "../mae-pipeline/Output/html/index.html")`'
#' output:
#'   html_document:
#'    code_folding: hide
#'    code_download: TRUE
#'---

saveRDS(snakemake, "tmp/outrider_overview.snakemake")
# snakemake <- readRDS("tmp/outrider_overview.snakemake")

summarie_titles <- paste(snakemake@config$GENE_ANNOTATION_NAMES, names(snakemake@config$outrider_filtered))
summaries <- paste('[', summarie_titles ,'](', gsub(snakemake@config$htmlOutputPath, ".", snakemake@input$summaries), ')', sep = '')
summaries <- paste(summaries, sep = '\n')
#' Summaries: `r summaries`
