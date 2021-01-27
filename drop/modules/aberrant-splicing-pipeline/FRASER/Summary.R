#'---
#' title: "FRASER Summary: `r gsub('_', ' ', snakemake@wildcards$dataset)`"
#' author: mumichae, vyepez, ischeller
#' wb:
#'  log:
#'    - snakemake: '`sm str(tmp_dir / "AS" / "{dataset}" / "FRASER_summary.Rds")`'
#'  params:
#'   - setup: '`sm cfg.AS.getWorkdir() + "/config.R"`'
#'   - workingDir: '`sm cfg.getProcessedDataDir() + "/aberrant_splicing/datasets/"`'
#'  input:
#'   - fdsin: '`sm cfg.getProcessedDataDir() + 
#'                 "/aberrant_splicing/datasets/savedObjects/{dataset}/" + 
#'                 "padjBetaBinomial_theta.h5"`'
#'   - results: '`sm cfg.getProcessedDataDir() + 
#'                   "/aberrant_splicing/results/{dataset}_results.tsv"`'
#'  output:
#'   - wBhtml: '`sm config["htmlOutputPath"] +
#'               "/AberrantSplicing/{dataset}_summary.html"`'
#'  type: noindex
#'---

saveRDS(snakemake, snakemake@log$snakemake)
source(snakemake@params$setup, echo=FALSE)

suppressPackageStartupMessages({
  library(cowplot)
})

#+ input
dataset    <- snakemake@wildcards$dataset
workingDir <- snakemake@params$workingDir

fds <- loadFraserDataSet(dir=workingDir, name=dataset)

#' Number of samples: `r nrow(colData(fds))`
#' 
#' Number of introns (psi5 or psi3): `r length(rowRanges(fds, type = "psi5"))`
#' 
#' Number of splice sites (theta): `r length(rowRanges(fds, type = "theta"))`

# used for most plots
dataset_title <- paste("Dataset:", dataset)


#' ## Hyperparameter optimization
for(type in psiTypes){
    g <- plotEncDimSearch(fds, type=type) 
    if (!is.null(g)) {
        g <- g + theme_cowplot(font_size = 16) + 
          ggtitle(paste0("Q estimation, ", type))
        print(g)
    }
}

#' ## Aberrantly spliced genes per sample
plotAberrantPerSample(fds, aggregate=TRUE, main=dataset_title) + 
    theme_cowplot(font_size = 16) +
    theme(legend.position = "top")

#' ## Batch Correlation: Samples x samples
topN <- 30000
topJ <- 10000
for(type in psiTypes){
    before <- plotCountCorHeatmap(
        fds,
        type = type,
        logit = TRUE,
        topN = topN,
        topJ = topJ,
        plotType = "sampleCorrelation",
        normalized = FALSE,
        annotation_col = NA,
        annotation_row = NA,
        sampleCluster = NA,
        plotMeanPsi=FALSE,
        plotCov = FALSE,
        annotation_legend = TRUE
    )
    before
    after <- plotCountCorHeatmap(
        fds,
        type = type,
        logit = TRUE,
        topN = topN,
        topJ = topJ,
        plotType = "sampleCorrelation",
        normalized = TRUE,
        annotation_col = NA,
        annotation_row = NA,
        sampleCluster = NA,
        plotMeanPsi=FALSE,
        plotCov = FALSE,
        annotation_legend = TRUE
    )
    after
}

#' # Results
res <- fread(snakemake@input$results)
file <- gsub(".html$", ".tsv", snakemake@output$wBhtml)
write_tsv(res, file=file)

#'
#' The results table can also be downloaded with the link below.
#+ echo=FALSE, results='asis'
cat(paste0("<a href='./", basename(file), "'>Download results table</a>"))

# round numbers
if(nrow(res) > 0){
  res[, pValue := signif(pValue, 3)]
  res[, padjust := signif(padjust, 3)]
  res[, deltaPsi := signif(deltaPsi, 2)]
  res[, zscore := signif(zScore, 2)]
  res[, psiValue := signif(psiValue, 2)]
  res[, pValueGene := signif(pValueGene, 2)]
  res[, padjustGene := signif(padjustGene, 2)]
}

#' ## Results table
DT::datatable(res, options=list(scrollX=TRUE), escape=FALSE, filter = 'top')

#' ## Samples table
DT::datatable(as.data.table(colData(fds)), options=list(scrollX=TRUE))
