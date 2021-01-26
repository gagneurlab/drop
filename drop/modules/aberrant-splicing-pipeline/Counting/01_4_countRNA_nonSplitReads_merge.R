#'---
#' title: Merge Nonsplit Counts
#' author: Luise Schuller
#' wb:
#'  log:
#'    - snakemake: '`sm str(tmp_dir / "AS" / "{dataset}" / "01_4_nonSplitReadsMerge.Rds")`'
#'  params:
#'   - setup: '`sm cfg.AS.getWorkdir() + "/config.R"`'
#'   - workingDir: '`sm cfg.getProcessedDataDir() + "/aberrant_splicing/datasets"`'
#'  threads: 20
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
source(snakemake@params$setup, echo=FALSE)

dataset    <- snakemake@wildcards$dataset
workingDir <- snakemake@params$workingDir
params <- snakemake@config$aberrantSplicing

register(MulticoreParam(snakemake@threads))
# Limit number of threads for DelayedArray operations
setAutoBPPARAM(MulticoreParam(snakemake@threads))

# Read FRASER object
fds <- loadFraserDataSet(dir=workingDir, name=paste0("raw-", dataset))

# Read splice site coordinates from RDS
splitCounts_gRanges <- readRDS(snakemake@input$gRangesNonSplitCounts)

# If samples are recounted, remove the merged ones
nonSplitCountsDir <- file.path(workingDir, "savedObjects", 
                            paste0("raw-", dataset), 'nonSplitCounts')
if(params$recount == TRUE & dir.exists(nonSplitCountsDir)){
  unlink(nonSplitCountsDir, recursive = TRUE)
}

# Get and merge nonSplitReads for all sample ids
nonSplitCounts <- getNonSplitReadCountsForAllSamples(fds=fds, 
                                                     splitCountRanges=splitCounts_gRanges, 
                                                     minAnchor=5, 
                                                     recount=FALSE, 
                                                     longRead=params$longRead)

message(date(), ":", dataset, " nonSplit counts done")
