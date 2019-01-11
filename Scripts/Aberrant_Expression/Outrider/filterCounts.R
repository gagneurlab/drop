#'---
#' title: Filter Counts for OUTRIDER
#' author: Michaela Mueller
#' wb:
#'  input:
#'   - counts: '`sm config["PROC_RESULTS"] + "/{annotation}/counts/total_counts.Rds"`'
#'   - txdb: '`sm config["PROC_RESULTS"] + "/{annotation}/txdb.Rds"`'
#'  output:
#'   - filtered_counts: '`sm config["PROC_RESULTS"] + "/{annotation}/counts/filtered_counts.Rds"`'
#'   - ods: '`sm config["PROC_RESULTS"] + "/{annotation}/outrider/ods_unfitted.Rds"`'
#'   - plot: '`sm config["PROC_RESULTS"] + "/{annotation}/outrider/filtered_hist.png"`'
#'  type: script
#'---

suppressPackageStartupMessages({
    library(OUTRIDER)
    library(SummarizedExperiment)
    library(ggplot2)
    library(data.table)
    library(dplyr)
})

saveRDS(snakemake, "tmp/filter_counts.snakemake")
counts <- readRDS(snakemake@input$counts)
ods <- OutriderDataSet(counts)
colData(ods)$sampleID <- colnames(ods)

# TODO: Add batches to colData for heatmap
# colData(ods)$batch <- as.character(NA)
# for(i in seq_along(batches)){
#     idxSampleBatch <- colnames(ods) %in% colnames(ss_counts_gene[[i]])
#     colData(ods)$batch[idxSampleBatch] <- batches[i]
# }

# filter not expressed genes
gencode_txdb <- loadDb(snakemake@input$txdb)
seqlevelsStyle(gencode_txdb) <- "UCSC"
gencode_txdb <- keepStandardChromosomes(gencode_txdb)

ods <- filterExpression(ods, gtfFile=gencode_txdb, filter=FALSE)
g <- plotFPKM(ods) + theme_bw(base_size = 14)
ggsave(snakemake@output$plot, g)

ods <- filterExpression(ods, gtfFile=gencode_txdb, filter=TRUE, fpkmCutoff=snakemake@config$fpkmCutoff)

saveRDS(counts(ods), snakemake@output$filtered_counts)
saveRDS(ods, snakemake@output$ods)



