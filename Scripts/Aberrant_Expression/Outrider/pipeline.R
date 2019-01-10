#'---
#' title: Filter Counts for OUTRIDER
#' author: Michaela Mueller
#' wb:
#'  input:
#'   - counts: '`sm config["PROC_RESULTS"] + "/{annotation}/counts/total_counts.Rds"`'
#'   - txdb: '`sm config["PROC_RESULTS"] + "/{annotation}/txdb.Rds"`'
#'   - unique_gene_names: "resources/gencode_{annotation}_unique_gene_name.tsv"
#'  output:
#'   - filtered_counts: '`sm config["PROC_RESULTS"] + "/{annotation}/outrider/filtered_counts.Rds"`'
#'   - ods: '`sm config["PROC_RESULTS"] + "/{annotation}/outrider/ods.Rds"`'
#'   - plot: '`sm config["PROC_RESULTS"] + "/{annotation}/outrider/filtered_hist.png"`'
#'  type: script
#'  threads: 30
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
ods <- filterExpression(ods, gtfFile=gencode_txdb, filter=FALSE)
g <- plotFPKM(ods) + theme_bw(base_size = 14)
ggsave(snakemake@output$plot, g)

ods <- filterExpression(ods, gtfFile=gencode_txdb, filter=TRUE, fpkmCutoff=snakemake@config$fpkmCutoff)
saveRDS(counts(ods), snakemake@output$filtered_counts)

# Add genes metainfo
genes_dt <- fread(snakemake@input$unique_gene_names)
rowData(ods)$geneID = row.names(ods)
# left join preserves order
rowData(ods) = left_join(as.data.table(rowData(ods)), genes_dt[,.(gene_id, gene_name, gene_type)], by = c("geneID" = "gene_id_unique"))
rownames(ods) = rowData(ods)$gene_name

# TODO: put into seperate script
# OUTRIDER pipeline
ods <- estimateSizeFactors(ods)
pars <- c(seq(5, min(c(40, ncol(ods), nrow(ods))), 2), 50, 70)
ods <- findEncodingDim(ods, lnorm = T, BPPARAM = MulticoreParam(snakemake@threads), params = pars)
# TODO: check encoding dimension plot
ods <- OUTRIDER(ods, BPPARAM = MulticoreParam(snakemake@threads))

# ods <- autoCorrect(ods, q = 60)  # Felix recommended, q = Ngenes / 4
# ods <- fit(ods)
# ods <- computePvalues(ods)
# ods <- computeZscores(ods)


# do it if you have time and a big memory 
#plotQQ(ods, global=TRUE)

saveRDS(ods, snakemake@output$ods)



