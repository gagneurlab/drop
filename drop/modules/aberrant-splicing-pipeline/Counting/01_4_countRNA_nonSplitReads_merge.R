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
#'                          "/aberrant_splicing/datasets/cache/raw-local-{dataset}/gRanges_NonSplitCounts.rds"`'
#'  output:
###   - countsSS: '`sm cfg.getProcessedDataDir() +
###                   "/aberrant_splicing/datasets/savedObjects/raw-local-{dataset}/rawCountsSS.h5"`'
#'   - done:     '`sm cfg.getProcessedDataDir() + 
#'                "/aberrant_splicing/datasets/savedObjects/raw-local-{dataset}/merge_theta.done"`'
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
fds <- loadFraserDataSet(dir=workingDir, name=paste0("raw-local-", dataset))

# Read splice site coordinates from RDS
splitCounts_gRanges <- readRDS(snakemake@input$gRangesNonSplitCounts)

# Remove the cache file of merged counts from previous runs
nonSplitCountsDir <- file.path(workingDir, "savedObjects", 
                            paste0("raw-local-", dataset), 'nonSplitCounts')
if(dir.exists(nonSplitCountsDir)){
  unlink(nonSplitCountsDir, recursive = TRUE)
}

# Get and merge nonSplitReads for all sample ids
nonSplitCounts <- getNonSplitReadCountsForAllSamples(fds=fds, 
                                                     splitCountRanges=splitCounts_gRanges, 
                                                     minAnchor=5, 
                                                     recount=FALSE, 
                                                     longRead=params$longRead)

# remove cache dir as its not needed later
if(dir.exists(nonSplitCountsDir)){
    unlink(nonSplitCountsDir, recursive = TRUE)
}

message(date(), ":", dataset, " nonSplit counts done")

file.create(snakemake@output$done)
