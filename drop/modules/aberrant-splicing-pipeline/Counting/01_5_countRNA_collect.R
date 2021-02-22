#'---
#' title: Collect all counts to FRASER Object
#' author: Luise Schuller
#' wb:
#'  log:
#'    - snakemake: '`sm str(tmp_dir / "AS" / "{dataset}" / "01_5_collect.Rds")`'
#'  input:
#'   - countsJ:  '`sm cfg.getProcessedDataDir() + 
#'                    "/aberrant_splicing/datasets/savedObjects/raw-{dataset}/rawCountsJ.h5"`'
#'   - countsSS: '`sm cfg.getProcessedDataDir() + 
#'                    "/aberrant_splicing/datasets/savedObjects/raw-{dataset}/rawCountsSS.h5"`'
#'   - gRangesSplitCounts: '`sm cfg.getProcessedDataDir() + 
#'                          "/aberrant_splicing/datasets/cache/raw-{dataset}/gRanges_splitCounts.rds"`'
#'   - spliceSites: '`sm cfg.getProcessedDataDir() + 
#'                   "/aberrant_splicing/datasets/cache/raw-{dataset}/spliceSites_splitCounts.rds"`'
#'  output:
#'   - counting_done: '`sm cfg.getProcessedDataDir() + 
#'                "/aberrant_splicing/datasets/savedObjects/raw-{dataset}/counting.done" `'
#'  type: script
#'---

saveRDS(snakemake, snakemake@log$snakemake)
suppressPackageStartupMessages(library(FRASER))

dataset      <- snakemake@wildcards$dataset
j_counts     <- snakemake@input$countsJ
theta_counts <- snakemake@input$countsSS

# Read FRASER object
fds <- loadFraserDataSet(file=j_counts)
splitCounts_gRanges <- readRDS(snakemake@input$gRangesSplitCounts)
spliceSiteCoords <- readRDS(snakemake@input$spliceSites)

# Get splitReads and nonSplitRead counts in order to store them in FRASER object
splitCounts_h5 <- HDF5Array::HDF5Array(j_counts, "rawCountsJ")
splitCounts_se <- SummarizedExperiment(
  colData = colData(fds),
  rowRanges = splitCounts_gRanges,
  assays = list(rawCountsJ=splitCounts_h5)
)


nonSplitCounts_h5 <- HDF5Array::HDF5Array(theta_counts, "rawCountsSS")
nonSplitCounts_se <- SummarizedExperiment(
  colData = colData(fds),
  rowRanges = spliceSiteCoords,
  assays = list(rawCountsSS=nonSplitCounts_h5)
)

# Add Counts to FRASER object
fds <- addCountsToFraserDataSet(fds=fds, 
    splitCounts=splitCounts_se,
    nonSplitCounts=nonSplitCounts_se)

# Save final FRASER object 
fds <- saveFraserDataSet(fds)

file.create(snakemake@output$counting_done)
