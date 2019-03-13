#'---
#' title: Filter Counts for OUTRIDER
#' author: Michaela Mueller
#' wb:
#'  input:
#'  - ods: '`sm config["PROC_RESULTS"] + "/{annotation}/outrider/{dataset}/ods_unfitted.Rds"`'
#'  output:
#'   - ods: '`sm config["PROC_RESULTS"] + "/{annotation}/outrider/{dataset}/ods.Rds"`'
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

# OUTRIDER pipeline
ods <- ods[mcols(ods)$passedFilter,] 
    
ods <- estimateSizeFactors(ods)
pars <- c(seq(5, min(c(40, ncol(ods), nrow(ods))), 2), 50, 70)

ods <- findEncodingDim(ods, lnorm = T, BPPARAM = MulticoreParam(snakemake@threads), params = pars)
# TODO: check encoding dimension plot

ods <- OUTRIDER(ods, BPPARAM = MulticoreParam(snakemake@threads))
message("outrider fitting finished")

# ods <- autoCorrect(ods, q = 60)  # Felix recommended, q = Ngenes / 4
# ods <- fit(ods)
# ods <- computePvalues(ods)
# ods <- computeZscores(ods)


# do it if you have time and a big memory 
# plotQQ(ods, global=TRUE)

row.names(ods) <- rowData(ods)$gene_name_unique

op <- snakemake@output$ods
op_date <- paste0(strsplit(op, "\\.")[[1]][1], "-", format(Sys.time(), "%Y%m%d") , ".Rds")
saveRDS(ods, op_date)

file.link(op, op_date)
# saveRDS(ods, snakemake@output$ods)



