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

fdsLocal <- loadFraserDataSet(dir=workingDir, name=paste0("raw-local-", dataset))
fdsMerge <- loadFraserDataSet(dir=workingDir, name=paste0("raw-", dataset))

has_external <- !(all(is.na(fdsMerge@colData$SPLICE_COUNTS_DIR)) || is.null(fdsMerge@colData$SPLICE_COUNTS_DIR))
if(has_external){
    fdsMerge@colData$isExternal <- as.factor(!is.na(fdsMerge@colData$SPLICE_COUNTS_DIR))
}else{
    fdsMerge@colData$isExternal <- as.factor(FALSE)
}
devNull <- saveFraserDataSet(fdsMerge,dir=workingDir, name=paste0("raw-", dataset))


#' ## Number of samples:   
#' Local: `r sum(!as.logical(fdsMerge@colData$isExternal))`  
#' External: `r sum(as.logical(fdsMerge@colData$isExternal))`  
#' 
#' **Using external counts**  
#' External counts introduce some complexity into the problem of counting junctions
#' because it is unknown whether or not a junction is not counted (because there are no reads)
#' compared to filtered and not present due to legal/personal sharing reasons. As a result,
#' after merging the local (counted from BAM files) counts and the external counts, only the junctions that are 
#' present in both remain. As a result it is likely that the number of junctions will decrease after merging.
#' 
#' 
#' ### Number of introns (psi5 or psi3) before and after merging:  
#' Local: `r length(rowRanges(fdsLocal, type = "psi5"))`  
#' Merged: `r length(rowRanges(fdsMerge, type = "psi5"))`  
#' 
#' ### Number of splice sites (theta) before and after merging: 
#' Local: `r length(rowRanges(fdsLocal, type = "theta"))`  
#' Merged: `r length(rowRanges(fdsMerge, type = "theta"))`  
#' 

#' ### Comparison of local and external counts  
if(has_external){
    externalCountIDs <- colData(fdsMerge)[as.logical(colData(fdsMerge)[,"isExternal"]),"sampleID"]
    localCountIDs <- colData(fdsMerge)[!as.logical(colData(fdsMerge)[,"isExternal"]),"sampleID"]

    cts <- K(fdsMerge,"psi5")
    ctsLocal<- cts[,localCountIDs,drop=FALSE]
    ctsExt<- cts[,externalCountIDs,drop=FALSE]

    rowMeanLocal <- rowMeans(ctsLocal)
    rowMeanExt <- rowMeans(ctsExt)

    dt <- data.table("Mean counts of local samples" = rowMeanLocal,
                     "Mean counts of external samples" = rowMeanExt)
                 
    ggplot(dt,aes(x = `Mean counts of local samples`, y= `Mean counts of external samples`)) +
       geom_hex() + theme_cowplot(font_size = 16) +
	   theme_bw() + scale_x_log10() + scale_y_log10() + 
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
