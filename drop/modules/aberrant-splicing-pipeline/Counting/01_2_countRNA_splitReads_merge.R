#'---
#' title: Merge Split Counts
#' author: Luise Schuller
#' wb:
#'  log:
#'    - snakemake: '`sm str(tmp_dir / "AS" / "{dataset}" / "01_2_splitReadsMerge.Rds")`'
#'  threads: 10
#'  input:
#'   - sample_counts: '`sm lambda w: cfg.AS.getSplitCountFiles(w.dataset)`'
#'   - fds_init:    '`sm cfg.getProcessedDataDir() + 
#'                   "/aberrant_splicing/datasets/savedObjects/raw-{dataset}/fds-init.done"`'
#'  output:
#'   - countsJ:  '`sm cfg.getProcessedDataDir() + 
#'                    "/aberrant_splicing/datasets/savedObjects/raw-{dataset}/rawCountsJ.h5"`'
#'   - gRangesSplitCounts: '`sm cfg.getProcessedDataDir() + 
#'                   "/aberrant_splicing/datasets/cache/raw-{dataset}/gRanges_splitCounts.rds"`'
#'   - gRangesNonSplitCounts: '`sm cfg.getProcessedDataDir() + 
#'                   "/aberrant_splicing/datasets/cache/raw-{dataset}/gRanges_NonSplitCounts.rds"`'
#'   - spliceSites: '`sm cfg.getProcessedDataDir() + 
#'                   "/aberrant_splicing/datasets/cache/raw-{dataset}/spliceSites_splitCounts.rds"`'
#'  type: script
#'---

saveRDS(snakemake, snakemake@log$snakemake)
suppressPackageStartupMessages(library(FRASER))

# input
dataset <- snakemake@wildcards$dataset
fdsIn   <- file.path(dirname(snakemake@input$fds_init), "fds-object.RDS")
params  <- snakemake@config$aberrantSplicing
minExpressionInOneSample <- params$minExpressionInOneSample
BPPARAM <- MulticoreParam(snakemake@threads)

# Set number of threads including for DelayedArray operations
register(BPPARAM)
DelayedArray::setAutoBPPARAM(BPPARAM)

# Force writing HDF5 files
options(FRASER.maxSamplesNoHDF5=-1)
options(FRASER.maxJunctionsNoHDF5=-1)

# Read FRASER object
fds <- loadFraserDataSet(file=fdsIn)

# enforce remerging if this script has to be run
splitCountsDir <- file.path(workingDir(fds), "savedObjects", name(fds), "splitCounts")
if(dir.exists(splitCountsDir)){
    unlink(splitCountsDir, recursive=TRUE)
}

# Get and merge splitReads for all sample ids
splitCounts <- getSplitReadCountsForAllSamples(fds=fds, outDir=splitCountsDir, recount=FALSE)

# Extract, annotate and save granges
splitCountRanges <- rowRanges(splitCounts)

# Annotate granges from the split counts
splitCountRanges <- FRASER:::annotateSpliceSite(splitCountRanges)
saveRDS(splitCountRanges, snakemake@output$gRangesSplitCounts)

# Create ranges for non split counts
# Subset by minExpression
maxCount <- rowMaxs(assay(splitCounts, "rawCountsJ"))
passed <- maxCount >= minExpressionInOneSample
# extract granges after filtering
splitCountRanges <- splitCountRanges[passed,]

saveRDS(splitCountRanges, snakemake@output$gRangesNonSplitCounts)

# Extract splitSiteCoodinates: extract donor and acceptor sites
# take filtered splicesite map
spliceSiteCoords <- FRASER:::extractSpliceSiteCoordinates(splitCountRanges, fds)
saveRDS(spliceSiteCoords, snakemake@output$spliceSites)


message(date(), ": ", dataset, " total no. splice junctions = ", 
        length(splitCounts))
