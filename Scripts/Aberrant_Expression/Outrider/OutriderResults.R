#'---
#' title: Analyze OUTRIDER Object
#' author: Michaela Mueller
#' wb:
#'  input:
#'   - ods: '`sm expand(config["PROC_RESULTS"] + "/{annotation}/outrider/ods.Rds", annotation=config["ANNOTATIONS"])`'
#' output: 
#'   html_document
#'---

suppressPackageStartupMessages({
    library(OUTRIDER)
    library(SummarizedExperiment)
    library(ggplot2)
    library(data.table)
    library(dplyr)
})

saveRDS(snakemake, "tmp/outrider_analysis.snakemake")
ods <- sapply(snakemake@input$ods, readRDS)

#' # Quality Control
barplot(sizeFactors(ods), main = "Size Factors")

#+ heatmap, fig.height=8, fig.width=8
plotCountCorHeatmap(ods)

res_ss <- results(ods_ss)
res_ss[, IS_RNA_SEQ_STRANDED := T]


ods_nss <- readRDS("/s/project/genetic_diagnosis/processed_results/ods_batches0_1_th_nss.Rds")
res_nss <- results(ods_nss)
res_nss[, IS_RNA_SEQ_STRANDED := F]

res <- rbind(res_ss, res_nss)

res[, LAB := "PROKISCH"]
res[sampleID %in%  sat[, ID_Links], LAB := "HAACK"]


# do some global plotting
ods <- plotCountCorHeatmap(ods, normalized=FALSE, rowCoFactor="batch")
# ods <- ods[,!colData(ods)$clusterNumber %in% c("3", "4")]

ods <- plotCountCorHeatmap(ods, normalized=TRUE, rowCoFactor="batch")

png("/s/project/genetic_diagnosis/processed_results/barplot_aberrant_genes_per_sample.png", width=900, height=700, res=120)
plotAberrantPerSample(ods)
dev.off()