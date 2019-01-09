#'---
#' title: Run OUTRIDER Pipeline
#' author: Michaela Mueller
#' wb:
#'  input:
#'   - counts: '`sm config["PROC_RESULTS"] + "/counts/{annotation}/total_counts.Rds"`'
#'   - txdb: '`sm config["PROC_RESULTS"] + "/counts/{annotation}/txdb.Rds"`'
#'   - unique_gene_names: "resources/gencode_{annotation}_unique_gene_name.tsv"
#'  output:
#'   - ods: '`sm config["PROC_RESULTS"] + "/counts/{annotation}/ods.Rds"`'
#'   - plot: '`sm config["PROC_RESULTS"] + "/counts/{annotation}/filtered_hist.png"`'
#'  type: script
#'---

suppressPackageStartupMessages({
    library(OUTRIDER)
    library(SummarizedExperiment)
    library(ggplot2)
    library(data.table)
    library(dplyr)
})

saveRDS(snakemake, "tmp/outrider.snakemake")
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

ods <- filterExpression(ods, gtfFile=gencode_txdb, filter=TRUE)

saveRDS(ods, snakemake@output$ods)


# Add genes metainfo
genes_dt <- fread(snakemake@input$unique_gene_names)
gene_name_as_row_name <- function(ods, genes_dt){
    rowData(ods)$geneID = row.names(ods)
    # left join preserves order
    rowData(ods) = left_join(as.data.table(rowData(ods)), genes_dt[,.(gene_id, gene_name, gene_type)], by = c("geneID" = "gene_id_unique"))
    rownames(ods) = rowData(ods)$gene_name
    ods
}

ods <- gene_name_as_row_name(ods)


# run full outrider
ods <- estimateSizeFactors(ods)
# saveRDS(counts(ods, normalized = T), "/s/project/genetic_diagnosis/processed_data/Rds/batches2_3_4_counts_ss_fpkmfiltered_normsizefactors.Rds")
saveRDS(ods, "~/Downloads/ods.Rds")
ods <- readRDS("~/Downloads/ods.Rds")

pars <- c(seq(5, min(c(40, ncol(ods), nrow(ods))), 2), 50, 70)
ods <- findEncodingDim(ods, lnorm = T, BPPARAM = MulticoreParam(30), params = pars)
# TODO: check encoding dimension plot
ods <- OUTRIDER(ods, BPPARAM = MulticoreParam(30))

# ods <- autoCorrect(ods, q = 60)  # Felix recommended, q = Ngenes / 4
# ods <- fit(ods)
# ods <- computePvalues(ods)
# ods <- computeZscores(ods)

library(biomaRt)
ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
ensgID <- rownames(ods)
ensg_match <- data.table(
    ensgIDv = ensgID,
    ensgID = gsub("\\.[0-9]+", "", ensgID)
)
ensg2match <- as.data.table(getBM(c("ensembl_gene_id", "hgnc_symbol"), filters = "ensembl_gene_id", values = ensg_match$ensgID, mart = ensembl))
ensg2match
ensg_match
res_merge <- merge(ensg2match, ensg_match, by.x="ensembl_gene_id", by.y="ensgID")
res_merge <- res_merge[order(ensembl_gene_id, hgnc_symbol)][!duplicated(ensgIDv)]
setkey(res_merge, "ensgIDv")
res_merge
toreplace <- res_merge[rownames(ods)]

rownames(ods) <- toreplace[,hgnc_symbol]

rownames(ods)[toreplace[,hgnc_symbol == "" | is.na(hgnc_symbol)]] <- toreplace[hgnc_symbol == "" | is.na(hgnc_symbol), ensgIDv]

# saveRDS(ods, "/s/project/genetic_diagnosis/processed_results/ods_4batches.Rds")
saveRDS(ods_ss, "/s/project/genetic_diagnosis/processed_results/ods_batches2_3_4_th_ss.Rds")
saveRDS(ods_nss, "/s/project/genetic_diagnosis/processed_results/ods_batches0_1_th_nss.Rds")

ods_ss <- readRDS("/s/project/genetic_diagnosis/processed_results/ods_batches2_3_4_th_ss.Rds")

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

# do it if you have time and a big memory 
# plotQQ(ods, global=TRUE)

