#'---
#' title: Merge Split Counts
#' author: Luise Schuller
#' wb:
#'  py:
#'  - |
#'   def getSplitCountFiles(dataset):
#'       ids = parser.fraser_ids[dataset]
#'       file_stump = parser.getProcDataDir() + f"/aberrant_splicing/datasets/cache/raw-{dataset}/sample_tmp/splitCounts/"
#'       return expand(file_stump + "sample_{sample_id}.done", sample_id=ids) 
#'  params:
#'   - tmpdir: '`sm drop.getMethodPath(METHOD, "tmp_dir")`'
#'   - workingDir: '`sm parser.getProcDataDir() + "/aberrant_splicing/datasets"`'
#'  threads: 20
#'  input:
#'   - sample_counts: '`sm lambda wildcards: getSplitCountFiles(wildcards.dataset)`'
#'  output:
#'   - countsJ: '`sm parser.getProcDataDir() +
#'                   "/aberrant_splicing/datasets/savedObjects/raw-{dataset}/rawCountsJ.h5"`'
#'   - gRangesSplitCounts: '`sm parser.getProcDataDir() + 
#'                          "/aberrant_splicing/datasets/cache/raw-{dataset}/gRanges_splitCounts.rds"`'
#'   - gRangesNonSplitCounts: '`sm parser.getProcDataDir() + 
#'                          "/aberrant_splicing/datasets/cache/raw-{dataset}/gRanges_NonSplitCounts.rds"`'
#'   - spliceSites: '`sm parser.getProcDataDir() + 
#'                   "/aberrant_splicing/datasets/cache/raw-{dataset}/spliceSites_splitCounts.rds"`'
#'  type: script
#'---
saveRDS(snakemake, file.path(snakemake@params$tmpdir, "FRASER_01_2.snakemake"))
# snakemake <- readRDS(".drop/tmp/AS/FRASER_01_2.snakemake")

source("Scripts/_helpers/config.R")

dataset    <- snakemake@wildcards$dataset
workingDir <- snakemake@params$workingDir
params <- snakemake@config$aberrantSplicing
minExpressionInOneSample <- params$minExpressionInOneSample


register(MulticoreParam(snakemake@threads))
# Limit number of threads for DelayedArray operations
setAutoBPPARAM(MulticoreParam(snakemake@threads))

# Read FRASER object
fds <- loadFraserDataSet(dir=workingDir, name=paste0("raw-", dataset))

# If samples are recounted, remove the merged ones
splitCountsDir <- file.path(workingDir, "savedObjects", 
                            paste0("raw-", dataset), 'splitCounts')
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

# Create ranges for non split counts
# Subset by minExpression
maxCount <- rowMaxs(assay(splitCounts, "rawCountsJ"))
passed <- maxCount >= minExpressionInOneSample
# extract granges after filtering
splitCountRanges <- splitCountRanges[passed,]

saveRDS(splitCountRanges, snakemake@output$gRangesNonSplitCounts)

# Extract splitSiteCoodinates: extract donor and acceptor sites
# take either filtered or full fds
spliceSiteCoords <- FRASER:::extractSpliceSiteCoordinates(splitCountRanges, fds)
saveRDS(spliceSiteCoords, snakemake@output$spliceSites)


message(date(), ": ", dataset, " total no. splice junctions = ", 
        length(splitCounts))
