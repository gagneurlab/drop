#'---
#' title: Filter Counts for OUTRIDER
#' author: Michaela Mueller
#' wb:
#'  params:
#'   - tmpdir: '`sm drop.getMethodPath(METHOD, "tmp_dir")`'
#'  input:
#'   - counts: '`sm parser.getProcDataDir() +
#'              "/aberrant_expression/{annotation}/outrider/{dataset}/total_counts.Rds"`'
#'   - txdb: '`sm parser.getProcDataDir() +
#'            "/aberrant_expression/{annotation}/txdb.db"`'
#'  output:
#'   - ods: '`sm parser.getProcResultsDir() +
#'           "/aberrant_expression/{annotation}/outrider/{dataset}/ods_unfitted.Rds"`'
#'  type: script
#'---

saveRDS(snakemake,  file.path(snakemake@params$tmpdir, "filter_counts.snakemake") )
# snakemake <- readRDS(".drop/tmp/AE/filter_counts.snakemake")

suppressPackageStartupMessages({
    library(data.table)
    library(GenomicFeatures)
    library(SummarizedExperiment)
    library(OUTRIDER)
})

counts <- readRDS(snakemake@input$counts)
colData(counts)$sampleID <- colData(counts)$RNA_ID
ods <- OutriderDataSet(counts)
txdb <- loadDb(snakemake@input$txdb)

# filter not expressed genes
fpkmCutoff <- snakemake@config$aberrantExpression$fpkmCutoff
ods <- filterExpression(ods, gtfFile=txdb, filter=FALSE,
                        fpkmCutoff=fpkmCutoff, addExpressedGenes=TRUE)

# add column for genes with at least 1 gene
rowData(ods)$counted1sample = rowSums(assay(ods)) > 0

# Save the ods before filtering to preserve the original number of genes
saveRDS(ods, snakemake@output$ods)
