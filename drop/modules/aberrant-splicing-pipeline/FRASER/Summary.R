#'---
#' title: "FRASER Summary: `r paste(snakemake@wildcards$dataset, snakemake@wildcards$annotation, sep = '--')`"
#' author: mumichae, vyepez, ischeller
#' wb:
#'  log:
#'    - snakemake: '`sm str(tmp_dir / "AS" / "{dataset}--{annotation}" / "FRASER_summary.Rds")`'
#'  params:
#'   - setup: '`sm cfg.AS.getWorkdir() + "/config.R"`'
#'  input:
#'   - fdsin: '`sm cfg.getProcessedResultsDir() + 
#'                 "/aberrant_splicing/datasets/merged/savedObjects/{dataset}--{annotation}/fds-object.RDS"`'
#'   - results: '`sm cfg.getProcessedResultsDir() + 
#'                   "/aberrant_splicing/results/{annotation}/fraser/{dataset}/results.tsv"`'
#'  output:
#'   - wBhtml: '`sm config["htmlOutputPath"] +
#'               "/AberrantSplicing/{dataset}--{annotation}_summary.html"`'
#'   - res_html: '`sm config["htmlOutputPath"] +
#'               "/AberrantSplicing/FRASER_results_{dataset}--{annotation}.tsv"`'
#'  type: noindex
#'---

saveRDS(snakemake, snakemake@log$snakemake)
source(snakemake@params$setup, echo=FALSE)

suppressPackageStartupMessages({
  library(cowplot)
})

#+ input
dataset    <- snakemake@wildcards$dataset
annotation <- snakemake@wildcards$annotation
padj_cutoff <- snakemake@config$aberrantSplicing$padjCutoff
zScore_cutoff <- snakemake@config$aberrantSplicing$zScoreCutoff
deltaPsi_cutoff <- snakemake@config$aberrantSplicing$deltaPsiCutoff


fds <- loadFraserDataSet(file=snakemake@input$fdsin)

#' Number of samples: `r nrow(colData(fds))`
#' 
#' Number of introns (psi5 or psi3): `r length(rowRanges(fds, type = "psi5"))`
#' 
#' Number of splice sites (theta): `r length(rowRanges(fds, type = "theta"))`

# used for most plots
dataset_title <- paste0("Dataset: ", dataset, "--", annotation)


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
plotAberrantPerSample(fds, padjCutoff = padj_cutoff, zScoreCutoff = zScore_cutoff, deltaPsiCutoff = deltaPsi_cutoff,
                      aggregate=TRUE, main=dataset_title) + 
  theme_cowplot(font_size = 16) +
  theme(legend.position = "top")

#' ## Batch Correlation: samples x samples
topN <- 30000
topJ <- 10000
for(type in psiTypes){
  before <- plotCountCorHeatmap(
    object=fds,
    type = type,
    logit = TRUE,
    topN = topN,
    topJ = topJ,
    plotType = "sampleCorrelation",
    normalized = FALSE,
    annotation_col = "isExternal",
    annotation_row = NA,
    sampleCluster = NA,
    minDeltaPsi = snakemake@config$aberrantSplicing$minDeltaPsi,
    plotMeanPsi=FALSE,
    plotCov = FALSE,
    annotation_legend = TRUE
  )
  before
  after <- plotCountCorHeatmap(
    object=fds,
    type = type,
    logit = TRUE,
    topN = topN,
    topJ = topJ,
    plotType = "sampleCorrelation",
    normalized = TRUE,
    annotation_col = "isExternal",
    annotation_row = NA,
    sampleCluster = NA,
    minDeltaPsi = snakemake@config$aberrantSplicing$minDeltaPsi,
    plotMeanPsi=FALSE,
    plotCov = FALSE,
    annotation_legend = TRUE
  )
  after
}

#' ## Results
res <- fread(snakemake@input$results)
file <- snakemake@output$res_html
write_tsv(res, file = file)
#+ echo=FALSE, results='asis'
cat(paste0("<a href='./", basename(file), "'>Download FRASER results table</a>"))

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

DT::datatable(
  head(res, 1000),
  caption = 'FRASER results (up to 1,000 rows shown)',
  options=list(scrollX=TRUE),
  escape=FALSE,
  filter = 'top'
)

