#'---
#' title: "Count Summary: `r gsub('_', ' ', snakemake@wildcards$dataset)`"
#' author: Christian Mertes
#' wb:
#'  log:
#'    - snakemake: '`sm str(tmp_dir / "AS" / "{dataset}" / "CountSummary.Rds")`'
#'  params:
#'   - setup: '`sm cfg.AS.getWorkdir() + "/config.R"`'
#'   - workingDir: '`sm cfg.getProcessedDataDir() + "/aberrant_splicing/datasets/"`'
#'  input:
#'   - filter: '`sm cfg.getProcessedDataDir() + 
#'                "/aberrant_splicing/datasets/savedObjects/{dataset}/filter.done" `'
#'  output:
#'   - wBhtml: '`sm config["htmlOutputPath"] + 
#'                  "/AberrantSplicing/{dataset}_countSummary.html"`'
#'  type: noindex
#'---

saveRDS(snakemake, snakemake@log$snakemake)
source(snakemake@params$setup, echo=FALSE)

suppressPackageStartupMessages({
  library(cowplot)
})
# opts_chunk$set(fig.width=10, fig.height=8)

#+ input
dataset    <- snakemake@wildcards$dataset
workingDir <- snakemake@params$workingDir

fds <- loadFraserDataSet(dir=workingDir, name=paste0("raw-", dataset))

#' Number of samples: `r nrow(colData(fds))`
#' 
#' Number of introns (psi5 or psi3): `r length(rowRanges(fds, type = "psi5"))`
#' 
#' Number of splice sites (theta): `r length(rowRanges(fds, type = "theta"))`
#' 
#' Introns that passed filter
table(mcols(fds, type="j")[,"passed"])

#' ## Expression filtering
#' Min expression cutoff: `r snakemake@config$aberrantSplicing$minExpressionInOneSample`
plotFilterExpression(fds) + theme_cowplot(font_size = 16)

#' ## Variability filtering
#' Variability cutoff: `r snakemake@config$aberrantSplicing$minDeltaPsi`
plotFilterVariability(fds) + theme_cowplot(font_size = 16)

