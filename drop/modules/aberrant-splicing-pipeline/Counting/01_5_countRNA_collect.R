#'---
#' title: Collect all counts to FRASER Object
#' author: Luise Schuller
#' wb:
#'  log:
#'    - snakemake: '`sm str(tmp_dir / "AS" / "{dataset}" / "01_5_collect.Rds")`'
#'  params:
#'   - setup: '`sm cfg.AS.getWorkdir() + "/config.R"`'
#'   - workingDir: '`sm cfg.getProcessedDataDir() + "/aberrant_splicing/datasets"`'
#'  input:
#'    - countsSSdone: '`sm cfg.getProcessedDataDir() + 
#'                          "/aberrant_splicing/datasets/savedObjects/raw-local-{dataset}/merge_theta.done"`'
#'    - gRangesSplitCounts: '`sm cfg.getProcessedDataDir() + 
#'                          "/aberrant_splicing/datasets/cache/raw-local-{dataset}/gRanges_splitCounts.rds"`'
#'    - spliceSites: '`sm cfg.getProcessedDataDir() + 
#'                          "/aberrant_splicing/datasets/cache/raw-local-{dataset}/spliceSites_splitCounts.rds"`'
#'  output:
#'   - counting_done: '`sm cfg.getProcessedDataDir() + 
#'                          "/aberrant_splicing/datasets/savedObjects/raw-local-{dataset}/counting.done" `'
#'  type: script
#'---

saveRDS(snakemake, snakemake@log$snakemake)
source(snakemake@params$setup, echo=FALSE)

dataset    <- snakemake@wildcards$dataset
workingDir <- snakemake@params$workingDir
saveDir    <- dirname(snakemake@input$countsSSdone)

# Read FRASER object
fds <- loadFraserDataSet(dir=workingDir, name=paste0("raw-local-", dataset))
splitCounts_gRanges <- readRDS(snakemake@input$gRangesSplitCounts)
spliceSiteCoords <- readRDS(snakemake@input$spliceSites)

# Get splitReads and nonSplitRead counts in order to store them in FRASER object
splitCounts_h5 <- HDF5Array::HDF5Array(file.path(saveDir, "rawCountsJ.h5"), "rawCountsJ")
splitCounts_se <- SummarizedExperiment(
  colData = colData(fds),
  rowRanges = splitCounts_gRanges,
  assays = list(rawCountsJ=splitCounts_h5)
)


nonSplitCounts_h5 <- HDF5Array::HDF5Array(file.path(saveDir, "rawCountsSS.h5"), "rawCountsSS")
nonSplitCounts_se <- SummarizedExperiment(
  colData = colData(fds),
  rowRanges = spliceSiteCoords,
  assays = list(rawCountsSS=nonSplitCounts_h5)
)

# Add Counts to FRASER object
fds <- addCountsToFraserDataSet(fds=fds, splitCounts=splitCounts_se,
                                nonSplitCounts=nonSplitCounts_se)

# Save final FRASER object 
fds <- saveFraserDataSet(fds)

file.create(snakemake@output$counting_done)
