#'---
#' title: Merge Nonsplit Counts
#' author: Luise Schuller
#' wb:
#'  log:
#'    - snakemake: '`sm str(tmp_dir / "AS" / "{dataset}" / "01_4_nonSplitReadsMerge.Rds")`'
#'  threads: 10
#'  input:
#'   - sample_counts:  '`sm lambda w: cfg.AS.getNonSplitCountFiles(w.dataset)`'
#'   - gRangesNonSplitCounts: '`sm cfg.getProcessedDataDir() + 
#'                          "/aberrant_splicing/datasets/cache/raw-{dataset}/gRanges_NonSplitCounts.rds"`'
#'  output:
#'   - countsSS: '`sm cfg.getProcessedDataDir() +
#'                   "/aberrant_splicing/datasets/savedObjects/raw-{dataset}/rawCountsSS.h5"`'
#'  type: script
#'---

saveRDS(snakemake, snakemake@log$snakemake)
suppressPackageStartupMessages(library(FRASER))

dataset <- snakemake@wildcards$dataset
fds_in  <- file.path(dirname(snakemake@output$countsSS), "fds-object.RDS")
params  <- snakemake@config$aberrantSplicing
BPPARAM <- MulticoreParam(snakemake@threads)

# Set number of threads including for DelayedArray operations
register(BPPARAM)
DelayedArray::setAutoBPPARAM(BPPARAM)

# Read FRASER object
fds <- loadFraserDataSet(file=fds_in)

# Read splice site coordinates from RDS
splitCounts_gRanges <- readRDS(snakemake@input$gRangesNonSplitCounts)

# If samples are recounted, remove the merged ones
nonSplitCountsDir <- file.path(workingDir(fds), "savedObjects", 
        name(fds), 'nonSplitCounts')
if(dir.exists(nonSplitCountsDir)){
  unlink(nonSplitCountsDir, recursive=TRUE)
}

# Get and merge nonSplitReads for all sample ids
nonSplitCounts <- getNonSplitReadCountsForAllSamples(fds=fds, 
        splitCountRanges=splitCounts_gRanges, 
        minAnchor=5, 
        recount=FALSE, 
        longRead=params$longRead)

message(date(), ":", dataset, " nonSplit counts done")
