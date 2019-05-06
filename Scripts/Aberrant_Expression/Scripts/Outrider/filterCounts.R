#'---
#' title: Filter Counts for OUTRIDER
#' author: Michaela Mueller
#' wb:
#'  input:
#'   - counts: '`sm config["PROC_RESULTS"] + "/{annotation}/counts/{dataset}/total_counts.Rds"`'
#'   - txdb: '`sm config["PROC_RESULTS"] + "/{annotation}/txdb.db"`'
#'  output:
#'   - ods: '`sm config["PROC_RESULTS"] + "/{annotation}/outrider/{dataset}/ods_unfitted.Rds"`'
#'   - plot: '`sm config["PROC_RESULTS"] + "/{annotation}/outrider/{dataset}/filtered_hist.png"`'
#'   - filtered_counts: '`sm config["PROC_RESULTS"] + "/{annotation}/counts/{dataset}/filtered_counts.Rds"`'
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
# snakemake <- readRDS("tmp/filter_counts.snakemake")
counts <- readRDS(snakemake@input$counts)
ods <- OutriderDataSet(counts)

# filter not expressed genes
gencode_txdb <- loadDb(snakemake@input$txdb)
seqlevelsStyle(gencode_txdb) <- "UCSC"
gencode_txdb <- keepStandardChromosomes(gencode_txdb)

ods <- filterExpression(ods, gtfFile=gencode_txdb, filter=FALSE, fpkmCutoff=snakemake@config$fpkmCutoff)
g <- plotFPKM(ods) + theme_bw(base_size = 14)
ggsave(snakemake@output$plot, g)

rowData(ods)$counted1sample = rowSums(assay(ods)) > 0

# Save the ods object before filtering, so as to preserve the original number of genes
saveRDS(ods, snakemake@output$ods)

# Save the filtered count matrix (as a matrix)
ods <- ods[mcols(ods)$passedFilter,]
saveRDS(counts(ods), snakemake@output$filtered_counts)



