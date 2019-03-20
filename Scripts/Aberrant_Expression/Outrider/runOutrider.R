#'---
#' title: Filter Counts for OUTRIDER
#' author: Michaela Mueller
#' wb:
#'  input:
#'   - ods: '`sm config["PROC_RESULTS"] + "/{annotation}/outrider/{dataset}/ods_unfitted.Rds"`'
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
    library(magrittr)
})

saveRDS(snakemake, "tmp/outrider.snakemake")
# snakemake <- readRDS("tmp/outrider.snakemake")
ods <- readRDS(snakemake@input$ods)

# OUTRIDER pipeline
ods <- ods[mcols(ods)$passedFilter,] 
    
ods <- estimateSizeFactors(ods)

a = 5 
b = min(ncol(ods), nrow(ods)) / 3   # N/3
Nsteps = min(20, ncol(ods)/3, nrow(ods)/3)   # Do either 20 steps or N
pars_q <- round(exp(seq(log(a),log(b),length.out = Nsteps))) %>% unique # Do unique in case 2 were repeated

ods <- findEncodingDim(ods, lnorm = T, BPPARAM = MulticoreParam(snakemake@threads), params = pars_q)

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

# Save the new ods with a date stamp
op_date <- paste0(strsplit(op, "\\.")[[1]][1], "-", format(Sys.time(), "%Y%m%d") , ".Rds")
saveRDS(ods, op_date)

# Create a link to the previous file
if(file.exists(op)) file.remove(op)
file.symlink(op_date, op)
