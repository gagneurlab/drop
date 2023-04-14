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
#'   - filter: '`sm expand(cfg.getProcessedDataDir() +
#'                  "/aberrant_splicing/datasets/savedObjects/{dataset}/filter_{version}.done", version=cfg.AS.get("FRASER_version"), allow_missing=True)`'
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
#' ## Number of introns:  
#' Local (before filtering): `r length(rowRanges(fdsLocal, type = "j"))`  
#' ```{asis, echo = has_external} 
#' Merged with external counts (before filtering):
#' ``` 
#' ```{r, eval = has_external, echo=FALSE} 
#' length(rowRanges(fdsMerge, type = "j"))
#' ```
#' After filtering: `r sum(mcols(fdsMerge, type="j")[,"passed"])`
#' 
#' ## Number of splice sites: 
#' Local: `r length(rowRanges(fdsLocal, type = "theta"))`  
#' ```{asis, echo = has_external}
#' Merged with external counts:
#' ``` 
#' ```{r, eval = has_external, echo=FALSE} 
#' length(rowRanges(fdsMerge, type = "theta"))
#' ```
#' 

#' ```{asis, echo = has_external}
#' ## Comparison of local and external counts  
#' **Using external counts**  
#' External counts introduce some complexity into the problem of counting junctions
#' because it is unknown whether or not a junction is not counted (because there are no reads)
#' compared to filtered and not present due to legal/personal sharing reasons. As a result,
#' after merging the local (counted from BAM files) counts and the external counts, only the junctions that are 
#' present in both remain. As a result it is likely that the number of junctions will decrease after merging.
#' 
#' ```
#' ```{r, eval = has_external, echo=has_external}
#' if(has_external){
#'     externalCountIDs <- colData(fdsMerge)[as.logical(colData(fdsMerge)[,"isExternal"]),"sampleID"]
#'     localCountIDs <- colData(fdsMerge)[!as.logical(colData(fdsMerge)[,"isExternal"]),"sampleID"]
#' 
#'     cts <- K(fdsMerge,"psi5")
#'     ctsLocal<- cts[,localCountIDs,drop=FALSE]
#'     ctsExt<- cts[,externalCountIDs,drop=FALSE]
#' 
#'     rowMeanLocal <- rowMeans(ctsLocal)
#'     rowMeanExt <- rowMeans(ctsExt)
#' 
#'     dt <- data.table("Mean counts of local samples" = rowMeanLocal,
#'                      "Mean counts of external samples" = rowMeanExt)
#'                  
#'     ggplot(dt,aes(x = `Mean counts of local samples`, y= `Mean counts of external samples`)) +
#'        geom_hex() + theme_cowplot(font_size = 16) +
#' 	   theme_bw() + scale_x_log10() + scale_y_log10() + 
#'        geom_abline(slope = 1, intercept =0) +
#'        scale_color_brewer(palette="Dark2") 
#' }
#' ```
#' 

#' ## Expression filtering
#' The expression filtering step removes introns that are lowly expressed. The requirements for an intron to pass this filter are:
#' 
#' * at least 1 sample has `r snakemake@config$aberrantSplicing$minExpressionInOneSample` counts (K) for the intron
#' * at least `r 100*(1-snakemake@config$aberrantSplicing$quantileForFiltering)`% of the samples need to have a total of at least `r snakemake@config$aberrantSplicing$quantileMinExpression` reads for the splice metric denominator (N) of the intron
plotFilterExpression(fdsMerge) + 
    labs(title="", x="Mean Intron Expression", y="Introns") +
    theme_cowplot(font_size = 16)

#' ## Variability filtering
#' The variability filtering step removes introns that have no or little variability in the splice metric values across samples. The requirement for an intron to pass this filter is:
#' 
#' * at least 1 sample has a difference of at least `r snakemake@config$aberrantSplicing$minDeltaPsi` in the splice metric compared to the mean splice metric of the intron
plotFilterVariability(fdsMerge) + 
    labs(title="", y="Introns") +
    theme_cowplot(font_size = 16)
