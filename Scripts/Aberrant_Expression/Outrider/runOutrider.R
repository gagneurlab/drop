#'---
#' title: Filter Counts for OUTRIDER
#' author: Michaela Mueller
#' wb:
#'  input:
#'  - ods: '`sm config["PROC_RESULTS"] + "/{annotation}/outrider/{strand}/ods_unfitted.Rds"`'
#'  output:
#'   - ods: '`sm config["PROC_RESULTS"] + "/{annotation}/outrider/{strand}/ods.Rds"`'
#'   - plot: '`sm config["PROC_RESULTS"] + "/{annotation}/outrider/{strand}/qqplot.png"`'
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

# ods <- autoCorrect(ods, q = 60)  # Felix recommended, q = Ngenes / 4
# ods <- fit(ods)
# ods <- computePvalues(ods)
# ods <- computeZscores(ods)


# do it if you have time and a big memory 
plotQQ(ods, global=TRUE)

row.names(ods) <- rowData(ods)$gene_name_unique

saveRDS(ods, snakemake@output$ods)



