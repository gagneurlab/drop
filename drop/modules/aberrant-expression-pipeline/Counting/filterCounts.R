#'---
#' title: Filter Counts for OUTRIDER
#' author: Michaela Mueller
#' wb:
#'  log:
#'   - snakemake: '`sm str(tmp_dir / "AE" / "{annotation}" / "{dataset}" / "filter.Rds")`'
#'  input:
#'   - counts: '`sm cfg.getProcessedDataDir() +
#'              "/aberrant_expression/{annotation}/outrider/{dataset}/total_counts.Rds"`'
#'   - txdb: '`sm cfg.getProcessedDataDir() +
#'            "/preprocess/{annotation}/txdb.db"`'
#'  output:
#'   - ods: '`sm cfg.getProcessedResultsDir() +
#'           "/aberrant_expression/{annotation}/outrider/{dataset}/ods_unfitted.Rds"`'
#'  type: script
#'---

saveRDS(snakemake, snakemake@log$snakemake)

suppressPackageStartupMessages({
    library(data.table)
    library(GenomicFeatures)
    library(SummarizedExperiment)
    library(OUTRIDER)
})

counts <- readRDS(snakemake@input$counts)
ods <- OutriderDataSet(counts)
txdb <- loadDb(snakemake@input$txdb)

# filter not expressed genes
fpkmCutoff <- snakemake@config$aberrantExpression$fpkmCutoff
ods <- filterExpression(ods, gtfFile=txdb, filter=FALSE,
                        fpkmCutoff=fpkmCutoff, addExpressedGenes=TRUE)

# add column for genes with at least 1 gene
rowData(ods)$counted1sample = rowSums(assay(ods)) > 0

# add gene ranges
gr <- unlist(endoapply(rowRanges(ods), range))
rowData(ods) <- cbind(rowData(ods), DataFrame(as.data.table(gr)))

# Save the ods before filtering to preserve the original number of genes
saveRDS(ods, snakemake@output$ods)
