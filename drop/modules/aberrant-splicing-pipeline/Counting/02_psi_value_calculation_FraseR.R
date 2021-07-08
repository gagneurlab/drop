#'---
#' title: Calculate PSI values
#' author: Christian Mertes
#' wb:
#'  log:
#'   - snakemake: '`sm str(tmp_dir / "AS" / "{dataset}" / "02_PSIcalc.Rds")`'
#'  params:
#'   - workingDir: '`sm cfg.getProcessedDataDir() + "/aberrant_splicing/datasets/"`'
#'  threads: 30
#'  input:
#'   - setup: '`sm cfg.AS.getWorkdir() + "/config.R"`'
#'   - counting_done: '`sm cfg.getProcessedDataDir() + 
#'                "/aberrant_splicing/datasets/savedObjects/raw-{dataset}/counting.done" `'
#'  output:
#'  - theta:     '`sm cfg.getProcessedDataDir() +
#'                    "/aberrant_splicing/datasets/savedObjects/raw-{dataset}/theta.h5"`'
#'  type: script
#'--- 

saveRDS(snakemake, snakemake@log$snakemake)
source(snakemake@input$setup, echo=FALSE)

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
