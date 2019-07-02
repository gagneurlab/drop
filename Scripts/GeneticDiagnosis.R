#'---
#' title: Results Overview
#' author: mumichae
#' wb:
#'  input:
#'  - abb_expr: '`sm aberrantExp( "submodules/aberrant-expression-pipeline/Output/html/index.html")`'
#'  - abb_splicing: '`sm aberrantSplicing( "submodules/aberrant-splicing-pipeline/Output/html/index.html")`'
#'  - mae: '`sm mae( "submodules/mae-pipeline/Output/html/index.html")`'
#' output:
#'   html_document:
#'    code_folding: hide
#'    code_download: TRUE
#'---


# #'  - variants: '`sm variants( "submodules/variant-annotation-pipeline/Output/html/index.html")`'

saveRDS(snakemake, "tmp/outrider_overview.snakemake")
# snakemake <- readRDS("tmp/outrider_overview.snakemake")

gene_annotation_names <- names(snakemake@config$GENE_ANNOTATION)
summarie_titles <- paste(gene_annotation_names, names(snakemake@config$outrider_filtered))
summaries <- paste('[', summarie_titles ,'](', gsub(snakemake@config$htmlOutputPath, ".", snakemake@input$summaries), ')', sep = '')
summaries <- paste(summaries, sep = '\n')
#' Summaries: `r summaries`
