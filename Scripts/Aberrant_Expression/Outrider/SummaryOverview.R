#'---
#' title: Results Overview
#' author: mumichae
#' wb:
#'  input:
#'  - summaries: '`sm expand("Output/html/AberrantExpression/Outrider/{annotation}/OutriderSummary_{dataset}.html", annotation=config["ANNOTATIONS"] , dataset=[*config["outrider"]])`'
#' output:
#'   html_document:
#'    code_folding: hide
#'    code_download: TRUE
#'---

#  - summaries: '`sm expand("Output/html/AberrantExpression/Outrider/{annotation}/OutriderSummary_{dataset}.html", annotation=config["ANNOTATIONS"] , dataset=[*config["outrider_filtered"]])`'

saveRDS(snakemake, "tmp/outrider_overview.snakemake")
# snakemake <- readRDS("tmp/outrider_overview.snakemake")

summarie_titles <- paste(snakemake@config$ANNOTATIONS, names(snakemake@config$outrider))
summaries <- paste('[', summarie_titles ,'](', gsub(snakemake@config$htmlOutputPath, ".", snakemake@input$summaries), ')', sep = '')
summaries <- paste(summaries, sep = '\n')
#' Summaries: `r summaries`
