#'---
#' title: Filter Counts for OUTRIDER
#' author: Michaela Mueller
#' wb:
#'  input:
#'   - ods: '`sm config["PROC_RESULTS"] + "/{annotation}/outrider/ods_unfitted.Rds"`'
#'   - unique_gene_names: "/s/project/genetic_diagnosis/resource/gencode_{annotation}_unique_gene_name.tsv"
#'  output:
#'   - ods: '`sm config["PROC_RESULTS"] + "/{annotation}/outrider/ods.Rds"`'
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

saveRDS(snakemake, "tmp/outrider.snakemake")
# snakemake <- readRDS("tmp/outrider.snakemake")
ods <- readRDS(snakemake@input$ods)

# Add genes metainfo
genes_dt <- fread(snakemake@input$unique_gene_names)
rowData(ods)$geneID = row.names(ods)
# left join preserves order
rowData(ods) = left_join(as.data.table(rowData(ods)), genes_dt[,.(gene_id, gene_name, gene_type, gene_id_unique)], by = c("geneID" = "gene_id_unique"))
rownames(ods) = rowData(ods)$gene_name

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



