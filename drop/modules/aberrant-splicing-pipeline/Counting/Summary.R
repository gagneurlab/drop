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
#' **Using external counts**  
#' External counts introduce some complexity into the problem of counting junctions
#' because it is ambiguous whether or not a junction is not counted (because there are no reads)
#' compared to filtered and not present due to legal/personal sharing reasons. As a result,
#' after merging the local (fromBam) counts and the external counts, only the junctions that are exactly
#' the same in both remain. As a result it is likely that the number of junctions will decrease after a merge.
#' 
#' 
#' ### Number of introns (psi5 or psi3) before and after merging:  
#' Local (fromBam): `r length(rowRanges(fdsLocal, type = "psi5"))`  
#' Merged : `r length(rowRanges(fdsMerge, type = "psi5"))`  
#' 
#' ### Number of splice sites (theta) before and after merging: 
#' Local (fromBam): `r length(rowRanges(fdsLocal, type = "theta"))`  
#' Merged: `r length(rowRanges(fdsMerge, type = "theta"))`  
#' 

#' ### Comparison of local and external counts  
if(has_external){
externalCountIDs <- colData(fdsMerge)[colData(fdsMerge)[,"isExternal"],"sampleID"]
localCountIDs <- colData(fdsMerge)[!colData(fdsMerge)[,"isExternal"],"sampleID"]

cts <- K(fdsMerge,"psi5")
ctsLocal<- cts[,localCountIDs]
ctsExt<- cts[,externalCountIDs]

rowlgmLocal <- rowMeans(log(ctsLocal + 1))
rowlgmExt <- rowMeans(log(ctsExt + 1))

dt <- data.table("Local log mean counts" = rowlgmLocal,
                 "External log mean counts" = rowlgmExt)
                 
ggplot(dt,aes(x = `Local log mean counts`, y= `External log mean counts`)) +
   geom_point() + theme_cowplot(font_size = 16) +
   geom_abline(slope = 1, intercept =0) +
   scale_color_brewer(palette="Dark2") 
}else{
	print("No external counts, comparison is ommitted")
}

#' ## Expression filtering
#' Min expression cutoff: `r snakemake@config$aberrantSplicing$minExpressionInOneSample`
plotFilterExpression(fdsMerge) + theme_cowplot(font_size = 16)

#' ## Variability filtering
#' Variability cutoff: `r snakemake@config$aberrantSplicing$minDeltaPsi`
plotFilterVariability(fdsMerge) + theme_cowplot(font_size = 16)

#' Introns that passed filter (after merging)
table(mcols(fdsMerge, type="j")[,"passed"])
