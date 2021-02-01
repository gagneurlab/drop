#'---
#' title: Filter Counts for OUTRIDER
#' author: Michaela Mueller
#' wb:
#'  log:
#'   - snakemake: '`sm str(tmp_dir / "AE" / "{annotation}" / "{dataset}" / "runOUTRIDER.Rds")`'
#'  input:
#'   - ods: '`sm cfg.getProcessedResultsDir() + 
#'           "/aberrant_expression/{annotation}/outrider/{dataset}/ods_unfitted.Rds"`'
#'  output:
#'   - ods: '`sm cfg.getProcessedResultsDir() + 
#'           "/aberrant_expression/{annotation}/outrider/{dataset}/ods.Rds"`'
#'  type: script
#'  threads: 30
#'---

#+ echo=F
saveRDS(snakemake, snakemake@log$snakemake)

suppressPackageStartupMessages({
    library(OUTRIDER)
    library(SummarizedExperiment)
    library(ggplot2)
    library(data.table)
    library(dplyr)
    library(magrittr)
    library(tools)
})

ods <- readRDS(snakemake@input$ods)
implementation <- snakemake@config$aberrantExpression$implementation
mp <- snakemake@config$aberrantExpression$maxTestedDimensionProportion
register(MulticoreParam(snakemake@threads))

## subset filtered and estimate
ods <- ods[mcols(ods)$passedFilter,] 
ods <- estimateSizeFactors(ods)

## find optimal encoding dimension
a <- 5 
b <- min(ncol(ods), nrow(ods)) / mp   # N/3

maxSteps <- 15
if(mp < 4){
    maxSteps <- 20
}

Nsteps <- min(maxSteps, b)   # Do at most 20 steps or N/3
# Do unique in case 2 were repeated
pars_q <- round(exp(seq(log(a),log(b),length.out = Nsteps))) %>% unique
ods <- findEncodingDim(ods, params = pars_q, implementation = implementation)

## fit OUTRIDER
ods <- OUTRIDER(ods, implementation = implementation)
message("outrider fitting finished")

# Save the new ods with a date stamp
op <- snakemake@output$ods
op_date <- paste0(file_path_sans_ext(op), "-", format(Sys.time(), "%Y%m%d") , ".Rds")
saveRDS(ods, op_date)

# Create a link to the previous file
if(file.exists(op)) file.remove(op)
file.symlink(op_date, op)
