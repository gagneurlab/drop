#'---
#' title: Collect all counts to FRASER Object
#' author: Luise Schuller
#' wb:
#'  params:
#'   - tmpdir: '`sm drop.getMethodPath(METHOD, "tmp_dir")`'
#'   - workingDir: '`sm parser.getProcDataDir() + "/aberrant_splicing/datasets"`'
#'  input:
#'   - countsJ:  '`sm parser.getProcDataDir() + 
#'                    "/aberrant_splicing/datasets/savedObjects/raw-{dataset}/rawCountsJ.h5"`'
#'   - countsSS: '`sm parser.getProcDataDir() + 
#'                    "/aberrant_splicing/datasets/savedObjects/raw-{dataset}/rawCountsSS.h5"`'
#'   - gRangesSplitCounts: '`sm parser.getProcDataDir() + 
#'                          "/aberrant_splicing/datasets/cache/raw-{dataset}/gRanges_splitCounts.rds"`'
#'   - spliceSites: '`sm parser.getProcDataDir() + 
#'                   "/aberrant_splicing/datasets/cache/raw-{dataset}/spliceSites_splitCounts.rds"`'
#'  output:
#'   - counting_done: '`sm parser.getProcDataDir() + 
#'                "/aberrant_splicing/datasets/savedObjects/raw-{dataset}/counting.done" `'
#'  type: script
#'---
saveRDS(snakemake, file.path(snakemake@params$tmpdir, "FRASER_01_5.snakemake"))
# snakemake <- readRDS(".drop/tmp/AS/FRASER_01_5.snakemake")


source("Scripts/_helpers/config.R")

dataset    <- snakemake@wildcards$dataset
workingDir <- snakemake@params$workingDir


# Read FRASER object
fds <- loadFraserDataSet(dir=workingDir, name=paste0("raw-", dataset))
splitCounts_gRanges <- readRDS(snakemake@input$gRangesSplitCounts)
spliceSiteCoords <- readRDS(snakemake@input$spliceSites)

# Get splitReads and nonSplitRead counts in order to store them in FRASER object
splitCounts_h5 <- HDF5Array::HDF5Array(snakemake@input$countsJ, "rawCountsJ")
splitCounts_se <- SummarizedExperiment(
  colData = colData(fds),
  rowRanges = splitCounts_gRanges,
  assays = list(rawCountsJ=splitCounts_h5)
)


nonSplitCounts_h5 <- HDF5Array::HDF5Array(snakemake@input$countsSS, "rawCountsSS")
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
