#'---
#' title: Analyze OUTRIDER Object
#' author: Michaela Mueller, vyepez
#' wb:
#'  input:
#'   - ods_ss: '`sm expand(config["PROC_RESULTS"] + "/{annotation}/outrider/ss/ods.Rds", annotation=config["ANNOTATIONS"])`'
#'   - ods_ns: '`sm expand(config["PROC_RESULTS"] + "/{annotation}/outrider/ns/ods.Rds", annotation=config["ANNOTATIONS"])`'
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
# snakemake <- readRDS("tmp/outrider_analysis.snakemake")

ods <- sapply(snakemake@input$ods, readRDS)
names(ods)
ods <- ods[[3]]

dim(ods)

#' # Quality Control
barplot(sort(sizeFactors(ods)), main = "Size Factors", xaxt = 'n')

#+ heatmap, fig.height=8, fig.width=8
plotCountCorHeatmap(ods)

res <- results(ods)
res[, IS_RNA_SEQ_STRANDED := T]


# ods_nss <- readRDS("/s/project/genetic_diagnosis/processed_results/ods_batches0_1_th_nss.Rds")
# res_nss <- results(ods_nss)
# res_nss[, IS_RNA_SEQ_STRANDED := F]

# res <- rbind(res_ss, res_nss)

# res[sampleID %in%  sat[, ID_Links], LAB := "HAACK"]


# do some global plotting, add Batch to rowCoFactor
# ods <- plotCountCorHeatmap(ods, normalized=FALSE, rowCoFactor="batch")
# ods <- ods[,!colData(ods)$clusterNumber %in% c("3", "4")]
# ods <- plotCountCorHeatmap(ods, normalized=TRUE, rowCoFactor="batch")

png("/s/project/genetic_diagnosis/processed_results/barplot_aberrant_genes_per_sample.png", width=900, height=700, res=120)
plotAberrantPerSample(ods)
dev.off()