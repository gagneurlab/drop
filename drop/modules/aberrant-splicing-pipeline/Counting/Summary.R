#'---
#' title: "Count Summary: `r gsub('_', ' ', snakemake@wildcards$dataset)`"
#' author: Christian Mertes
#' wb:
#'  log:
#'    - snakemake: '`sm str(tmp_dir / "AS" / "{dataset}" / "CountSummary.Rds")`'
#'  params:
#'   - setup: '`sm cfg.AS.getWorkdir() + "/config.R"`'
#'   - workingDir: '`sm cfg.getProcessedDataDir() + "/aberrant_splicing/datasets/merged/"`'
#'   - workingDirLocal: '`sm cfg.getProcessedDataDir() + "/aberrant_splicing/datasets/fromBam/"`'
#'  input:
#'   - filter: '`sm cfg.getProcessedDataDir() + 
#'                "/aberrant_splicing/datasets/merged/savedObjects/{dataset}/filter.done" `'
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
workingDirLocal <- snakemake@params$workingDirLocal

fdsLocal <- loadFraserDataSet(dir=workingDirLocal, name=paste0("raw-", dataset))
fdsMerge <- loadFraserDataSet(dir=workingDir, name=paste0("raw-", dataset))

has_external <- !(all(is.na(fdsMerge@colData$SPLICE_COUNTS_DIR)) || is.null(fdsMerge@colData$SPLICE_COUNTS_DIR))
if(has_external){
    fdsMerge@colData$isExternal <- !is.na(fdsMerge@colData$SPLICE_COUNTS_DIR)
}else{
    fdsMerge@colData$isExternal <- FALSE
}

#' ## Number of samples:   
#' Local (fromBam): `r sum(!fdsMerge@colData$isExternal)`  
#' External: `r sum(fdsMerge@colData$isExternal)`  
#' 
#' ### Number of introns (psi5 or psi3):  
#' Local (fromBam): `r length(rowRanges(fdsLocal, type = "psi5"))`  
#' Merged : `r length(rowRanges(fdsMerge, type = "psi5"))`  
#' 
#' ### Number of splice sites (theta): 
#' Local (fromBam): `r length(rowRanges(fdsLocal, type = "theta"))`  
#' Merged: `r length(rowRanges(fdsMerge, type = "theta"))`  
#' 
#' Introns that passed filter (after merging)
table(mcols(fdsMerge, type="j")[,"passed"])

#' ## Expression filtering
#' Min expression cutoff: `r snakemake@config$aberrantSplicing$minExpressionInOneSample`
plotFilterExpression(fdsMerge) + theme_cowplot(font_size = 16)

#' ## Variability filtering
#' Variability cutoff: `r snakemake@config$aberrantSplicing$minDeltaPsi`
plotFilterVariability(fdsMerge) + theme_cowplot(font_size = 16)

