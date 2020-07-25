#'---
#' title: Nonsplit Counts
#' author: Luise Schuller
#' wb:
#'  params:
#'   - tmpdir: '`sm drop.getMethodPath(METHOD, "tmp_dir")`'
#'   - workingDir: '`sm parser.getProcDataDir() + "/aberrant_splicing/datasets"`'
#'  input:
#'   - spliceSites: '`sm parser.getProcDataDir() + 
#'                   "/aberrant_splicing/datasets/cache/raw-{dataset}/spliceSites_splitCounts.rds"`'
#'  output:
#'   - done_sample_nonSplitCounts : '`sm parser.getProcDataDir() + 
#'                   "/aberrant_splicing/datasets/cache/raw-{dataset}/sample_tmp/nonSplitCounts/sample_{sample_id}.done"`' 
#'  threads: 3
#'  type: script
#'---
saveRDS(snakemake, file.path(snakemake@params$tmpdir, "FRASER_01_3.snakemake"))
# snakemake <- readRDS(".drop/tmp/AS/FRASER_01_3.snakemake")

source("Scripts/_helpers/config.R")

dataset    <- snakemake@wildcards$dataset
colDataFile <- snakemake@input$colData
workingDir <- snakemake@params$workingDir
params <- snakemake@config$aberrantSplicing

# Read FRASER object
fds <- loadFraserDataSet(dir=workingDir, name=paste0("raw-", dataset))

# Get sample id from wildcard
sample_id <- snakemake@wildcards[["sample_id"]]


# Read splice site coordinates from RDS
spliceSiteCoords <- readRDS(snakemake@input$spliceSites)

# Count nonSplitReads for given sample id
sample_result <- countNonSplicedReads(sample_id,
                                      splitCountRanges = NULL,
                                      fds = fds,
                                      NcpuPerSample = snakemake@threads,
                                      minAnchor=5,
                                      recount=params$recount,
                                      spliceSiteCoords=spliceSiteCoords,
                                      longRead=params$longRead)

message(date(), ": ", dataset, ", ", sample_id,
        " no. splice junctions (non split counts) = ", length(sample_result))

file.create(snakemake@output$done_sample_nonSplitCounts)
