#'---
#' title: "OUTRIDER Summary: `r paste(snakemake@wildcards$dataset, snakemake@wildcards$annotation, sep = '--')`"
#' author: mumichae, vyepez
#' wb:
#'  log:
#'   - snakemake: '`sm str(tmp_dir / "AE" / "{annotation}" / "{dataset}" / "OUTRIDER_summary.Rds")`'
#'  params:
#'   - padjCutoff: '`sm cfg.AE.get("padjCutoff")`'
#'   - zScoreCutoff: '`sm cfg.AE.get("zScoreCutoff")`'
#'  input:
#'   - ods: '`sm cfg.getProcessedResultsDir() +
#'           "/aberrant_expression/{annotation}/outrider/{dataset}/ods.Rds"`'
#'   - results: '`sm cfg.getProcessedResultsDir() +
#'               "/aberrant_expression/{annotation}/outrider/{dataset}/OUTRIDER_results.tsv"`'
#'  output:
#'   - wBhtml: '`sm config["htmlOutputPath"] + 
#'              "/AberrantExpression/Outrider/{annotation}/Summary_{dataset}.html"`'
#'   - res_html: '`sm config["htmlOutputPath"] + 
#'              "/AberrantExpression/Outrider/{annotation}/OUTRIDER_results_{dataset}.tsv"`'
#'  type: noindex
#' output:
#'  html_document:
#'   code_folding: hide
#'   code_download: TRUE
#'---

#+ echo=F
saveRDS(snakemake, snakemake@log$snakemake)

suppressPackageStartupMessages({
  library(OUTRIDER)
  library(SummarizedExperiment)
  library(ggplot2)
  library(cowplot)
  library(data.table)
  library(dplyr)
  library(ggthemes)
})

# used for most plots
dataset_title <- paste("Dataset:", paste(snakemake@wildcards$dataset, snakemake@wildcards$annotation, sep = '--'))

ods <- readRDS(snakemake@input$ods)
if(is.null(colData(ods)$isExternal)) colData(ods)$isExternal <- FALSE

#' Number of samples: `r ncol(ods)`  
#' Number of expressed genes: `r nrow(ods)`  

#'
#' ## Visualize
#' ### Encoding dimension
plotEncDimSearch(ods) +
  labs(title = dataset_title) +
  theme_cowplot() +
  background_grid() +
  scale_color_brewer(palette="Dark2")


#' ### Aberrantly expressed genes per sample
plotAberrantPerSample(ods, main = dataset_title, 
                      padjCutoff = snakemake@params$padjCutoff,
                      zScoreCutoff = snakemake@params$zScoreCutoff)


#' ### Batch correction
#+ countCorHeatmap, fig.height=8, fig.width=8
plotCountCorHeatmap(ods, normalized = FALSE, colGroups = "isExternal", colColSet = "Dark2",
                    main = paste0('Raw Counts (', dataset_title, ')'))
plotCountCorHeatmap(ods, normalized = TRUE, colGroups = "isExternal", colColSet = "Dark2",
                    main = paste0('Normalized Counts (', dataset_title, ')'))


#' ### Expression by gene per sample
#+ geneSampleHeatmap, fig.height=12, fig.width=8
plotCountGeneSampleHeatmap(ods, normalized = FALSE, nGenes = 50, colGroups = "isExternal", colColSet = "Dark2",
                           main = paste0('Raw Counts (', dataset_title, ')'),
                           bcvQuantile = .95, show_names = 'row')
plotCountGeneSampleHeatmap(ods, normalized = TRUE, nGenes = 50, colGroups = "isExternal", colColSet = "Dark2",
                           main = paste0('Normalized Counts (',dataset_title,')'),
                           bcvQuantile = .95, show_names = 'row')


#' ### BCV - Biological coefficient of variation
# function to calculate BCV before autoencoder
estimateThetaWithoutAutoCorrect <- function(ods){
  
  ods1 <- OutriderDataSet(countData=counts(ods), colData=colData(ods))
  # use rowMeans as expected means
  normalizationFactors(ods1) <- matrix(rowMeans(counts(ods1)), 
                                       ncol=ncol(ods1), nrow=nrow(ods1))
  ods1 <- fit(ods1)
  theta(ods1)
  
  return(theta(ods1))
}

before <- data.table(when = "Before",
                     BCV = 1/sqrt(estimateThetaWithoutAutoCorrect(ods)))
after <- data.table(when = "After", BCV = 1/sqrt( theta( ods )))
bcv_dt <- rbind(before, after)

# boxplot of BCV Before and After Autoencoder
#+ BCV, fig.height=5, fig.width=6
ggplot(bcv_dt, aes(when, BCV)) +
  geom_boxplot() +
  theme_bw(base_size = 14) +
  labs(x = "Autoencoder correction", y = "Biological coefficient \nof variation",
       title = dataset_title)


#' ## Results
res <- fread(snakemake@input$results)

#' Total number of expression outliers: `r nrow(res)`  
#' Samples with at least one outlier gene: `r res[, uniqueN(sampleID)]`  

#'
#' ### Aberrant samples
#' 
#' An aberrant sample is one that has more than 0.1% of the total genes called as outliers.
if (nrow(res) > 0) {
  ab_table <- res[AberrantBySample > nrow(ods)/1000, .("Outlier genes" = .N), by = .(sampleID)] %>% unique
  if (nrow(ab_table) > 0) {
    setorder(ab_table, "Outlier genes") 
    DT::datatable(ab_table)
  } else {
    print("no aberrant samples")
  }
} else print("no aberrant samples")


#' ## Results table

## Save results table in the html folder and provide link to download
file <- snakemake@output$res_html
fwrite(res, file, sep = '\t', quote = F)


if(nrow(res) > 0){
  res[, pValue := format(pValue, scientific = T, digits = 3)]
  res[, padjust := format(padjust, scientific = T, digits = 3)]
  
  DT::datatable(head(res, 1000), caption = 'OUTRIDER results (up to 1,000 rows shown)',
                options=list(scrollX=TRUE), filter = 'top')
  
} else print("no significant results")

#+ echo=FALSE, results='asis'
cat(paste0("<a href='./", basename(file), "'>Download OUTRIDER results table</a>"))
