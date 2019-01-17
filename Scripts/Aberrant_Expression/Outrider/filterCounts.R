#'---
#' title: Filter Counts for OUTRIDER
#' author: Michaela Mueller
#' wb:
#'  input:
#'   - counts: '`sm config["PROC_RESULTS"] + "/{annotation}/counts/total_counts_{strand}.Rds"`'
#'   - txdb: '`sm config["PROC_RESULTS"] + "/{annotation}/txdb.db"`'
#'   - sample_anno: '`sm config["SAMPLE_ANNOTATION"]`'
#'  output:
#'   - ods: '`sm config["PROC_RESULTS"] + "/{annotation}/outrider/{strand}/ods_unfitted.Rds"`'
#'   - plot: '`sm config["PROC_RESULTS"] + "/{annotation}/outrider/{strand}/filtered_hist.png"`'
#'   - filtered_counts: '`sm config["PROC_RESULTS"] + "/{annotation}/counts/filtered_counts_{strand}.Rds"`'
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
# counts <- readRDS("/s/project/genetic_diagnosis/processed_results/v29/counts/total_counts.Rds")
ods <- OutriderDataSet(counts)
colData(ods)$sampleID <- colnames(ods)

# Add batches to colData for heatmap
sample_anno <- fread(snakemake@input$sample_anno)
cd <- colData(ods) %>% as.data.table
colData(ods)$batch <- left_join(cd, sample_anno[,.(RNA_ID, BATCH)], by = c("sampleID" = "RNA_ID"))$BATCH

# filter not expressed genes
gencode_txdb <- loadDb(snakemake@input$txdb)
seqlevelsStyle(gencode_txdb) <- "UCSC"
gencode_txdb <- keepStandardChromosomes(gencode_txdb)

ods <- filterExpression(ods, gtfFile=gencode_txdb, filter=FALSE)
g <- plotFPKM(ods) + theme_bw(base_size = 14)
ggsave(snakemake@output$plot, g)

rowData(ods)$counted1sample = rowSums(assay(ods)) > 0

# Save the ods object before filtering, so as to preserve the original number of genes
saveRDS(ods, snakemake@output$ods)

# Filter the ods object again and save the counts
ods <- filterExpression(ods, gtfFile=gencode_txdb, filter=TRUE, fpkmCutoff=snakemake@config$fpkmCutoff)
# Save the filtered count matrix (as a matrix)
saveRDS(counts(ods), snakemake@output$filtered_counts)



