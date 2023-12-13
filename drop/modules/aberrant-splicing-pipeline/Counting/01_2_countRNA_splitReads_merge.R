#'---
#' title: Merge Split Counts
#' author: Luise Schuller
#' wb:
#'  log:
#'    - snakemake: '`sm str(tmp_dir / "AS" / "{dataset}" / "01_2_splitReadsMerge.Rds")`'
#'  params:
#'   - setup: '`sm cfg.AS.getWorkdir() + "/config.R"`'
#'   - workingDir: '`sm cfg.getProcessedDataDir() + "/aberrant_splicing/datasets"`'
#'  threads: 20
#'  input:
#'   - sample_counts: '`sm lambda w: cfg.AS.getSplitCountFiles(w.dataset)`'
#'  output:
#'   - countsJ: '`sm cfg.getProcessedDataDir() +
#'                          "/aberrant_splicing/datasets/savedObjects/raw-local-{dataset}/rawCountsJ.h5"`'
#'   - gRangesSplitCounts: '`sm cfg.getProcessedDataDir() + 
#'                          "/aberrant_splicing/datasets/cache/raw-local-{dataset}/gRanges_splitCounts.rds"`'
#'   - gRangesNonSplitCounts: '`sm cfg.getProcessedDataDir() + 
#'                          "/aberrant_splicing/datasets/cache/raw-local-{dataset}/gRanges_NonSplitCounts.rds"`'
#'   - spliceSites: '`sm cfg.getProcessedDataDir() + 
#'                   "/aberrant_splicing/datasets/cache/raw-local-{dataset}/spliceSites_splitCounts.rds"`'
#'  type: script
#'---

saveRDS(snakemake, snakemake@log$snakemake)
source(snakemake@params$setup, echo=FALSE)

dataset    <- snakemake@wildcards$dataset
workingDir <- snakemake@params$workingDir
params <- snakemake@config$aberrantSplicing
minExpressionInOneSample <- params$minExpressionInOneSample


register(MulticoreParam(snakemake@threads))
# Limit number of threads for DelayedArray operations
setAutoBPPARAM(MulticoreParam(snakemake@threads))

# Read FRASER object
fds <- loadFraserDataSet(dir=workingDir, name=paste0("raw-local-", dataset))

# If samples are recounted, remove the merged ones
splitCountsDir <- file.path(workingDir, "savedObjects", 
                            paste0("raw-local-", dataset), 'splitCounts')
if(params$recount == TRUE & dir.exists(splitCountsDir)){
  unlink(splitCountsDir, recursive = TRUE)
}

# Get and merge splitReads for all sample ids
splitCounts <- getSplitReadCountsForAllSamples(fds=fds,
                                               recount=FALSE)
# Extract, annotate and save granges
splitCountRanges <- rowRanges(splitCounts)

# Annotate granges from the split counts
splitCountRanges <- FRASER:::annotateSpliceSite(splitCountRanges)
saveRDS(splitCountRanges, snakemake@output$gRangesSplitCounts)
# additionally save as tsv.gz (for easier AbSplice input)
fwrite(as.data.table(splitCountRanges), 
        gsub(".Rds", ".tsv.gz", snakemake@output$gRangesSplitCounts,
            ignore.case=TRUE))

# Create ranges for non split counts
# Subset by minExpression
maxCount <- rowMaxs(assay(splitCounts, "rawCountsJ"))
passed <- maxCount >= minExpressionInOneSample
# extract granges after filtering
splitCountRanges <- splitCountRanges[passed,]

saveRDS(splitCountRanges, snakemake@output$gRangesNonSplitCounts)

# Extract splitSiteCoodinates: extract donor and acceptor sites
# take either filtered or full fds
spliceSiteCoords <- FRASER:::extractSpliceSiteCoordinates(splitCountRanges)
saveRDS(spliceSiteCoords, snakemake@output$spliceSites)


message(date(), ": ", dataset, " total no. splice junctions = ", 
        length(splitCounts))
