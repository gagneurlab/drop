#'---
#' title: Calculate PSI values
#' author: Christian Mertes
#' wb:
#'  params:
#'   - tmpdir: '`sm drop.getMethodPath(METHOD, "tmp_dir")`'
#'   - workingDir: '`sm parser.getProcDataDir() + "/aberrant_splicing/datasets/"`'
#'  threads: 30
#'  input:
#'   - counting_done: '`sm parser.getProcDataDir() + 
#'                "/aberrant_splicing/datasets/savedObjects/raw-{dataset}/counting.done" `'
#'  output:
#'  - psiSS:     '`sm parser.getProcDataDir() + 
#'                    "/aberrant_splicing/datasets/savedObjects/raw-{dataset}/psiSite.h5"`'
#'  type: script
#'--- 

saveRDS(snakemake, file.path(snakemake@params$tmpdir, "FRASER_02.snakemake"))
# snakemake <- readRDS(".drop/tmp/AS/FRASER_02.snakemake")
source("Scripts/_helpers/config.R")

dataset    <- snakemake@wildcards$dataset
workingDir <- snakemake@params$workingDir
params <- snakemake@config$aberrantSplicing

register(MulticoreParam(snakemake@threads))
# Limit number of threads for DelayedArray operations
setAutoBPPARAM(MulticoreParam(snakemake@threads))

fds <- loadFraserDataSet(dir=workingDir, name=paste0("raw-", dataset))

# Calculating PSI values
fds <- calculatePSIValues(fds)

# FRASER object after PSI value calculation
fds <- saveFraserDataSet(fds)
