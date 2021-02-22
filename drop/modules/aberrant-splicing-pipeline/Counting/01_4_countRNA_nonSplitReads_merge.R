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

dataset  <- snakemake@wildcards$dataset
countsSS <- snakemake@output$countsSS
params <- snakemake@config$aberrantSplicing

register(MulticoreParam(snakemake@threads))
# Limit number of threads for DelayedArray operations
setAutoBPPARAM(MulticoreParam(snakemake@threads))

# Read FRASER object
fds <- loadFraserDataSet(file=file.path(dirname(countsSS), "fds-object.RDS"))

# Read splice site coordinates from RDS
splitCounts_gRanges <- readRDS(snakemake@input$gRangesNonSplitCounts)

# If samples are recounted, remove the merged ones
nonSplitCountsDir <- file.path(workingDir(fds), "savedObjects", 
        paste0("raw-", dataset), 'nonSplitCounts')
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
