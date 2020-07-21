#'---
#' title: Filter and clean dataset
#' author: Christian Mertes
#' wb:
#'  params:
#'   - tmpdir: '`sm drop.getMethodPath(METHOD, "tmp_dir")`'
#'   - workingDir: '`sm parser.getProcDataDir() + "/aberrant_splicing/datasets/"`'
#'  input:
#'   - psiSS:  '`sm parser.getProcDataDir()+ 
#'                  "/aberrant_splicing/datasets/savedObjects/raw-{dataset}/psiSite.h5"`'
#'  output:
#'   - fds: '`sm parser.getProcDataDir() +
#'                "/aberrant_splicing/datasets/savedObjects/{dataset}/fds-object.RDS"`'
#'   - done: '`sm parser.getProcDataDir() + 
#'                "/aberrant_splicing/datasets/savedObjects/{dataset}/filter.done" `'
#'  threads: 3
#'  type: script
#'---

saveRDS(snakemake, file.path(snakemake@params$tmpdir, "FRASER_03.snakemake"))
# snakemake <- readRDS(".drop/tmp/AS/FRASER_03.snakemake")

source("Scripts/_helpers/config.R")
opts_chunk$set(fig.width=12, fig.height=8)

# input
dataset    <- snakemake@wildcards$dataset
workingDir <- snakemake@params$workingDir
params <- snakemake@config$aberrantSplicing

fds <- loadFraserDataSet(dir=workingDir, name=paste0("raw-", dataset))

register(MulticoreParam(snakemake@threads))
# Limit number of threads for DelayedArray operations
setAutoBPPARAM(MulticoreParam(snakemake@threads))

# Apply filter
minExpressionInOneSample <- params$minExpressionInOneSample
minDeltaPsi <- params$minDeltaPsi

fds <- filterExpressionAndVariability(fds, 
                        minExpressionInOneSample = minExpressionInOneSample,
                        minDeltaPsi = minDeltaPsi,
                        filter=FALSE)
devNull <- saveFraserDataSet(fds)

# Keep junctions that pass filter
name(fds) <- dataset
if (params$filter == TRUE) {
    filtered <- mcols(fds, type="j")[,"passed"]
    fds <- fds[filtered,]
    message(paste("filtered to", nrow(fds), "junctions"))
}

fds <- saveFraserDataSet(fds)
file.create(snakemake@output$done)
