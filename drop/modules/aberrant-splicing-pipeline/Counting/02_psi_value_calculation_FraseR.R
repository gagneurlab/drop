#'---
#' title: Calculate PSI values
#' author: Christian Mertes
#' wb:
#'  log:
#'   - snakemake: '`sm str(tmp_dir / "AS" / "{dataset}" / "02_PSIcalc.Rds")`'
#'  threads: 10
#'  input:
#'   - counting_done: '`sm cfg.getProcessedDataDir() + 
#'                "/aberrant_splicing/datasets/savedObjects/raw-{dataset}/counting.done" `'
#'  output:
#'  - theta:     '`sm cfg.getProcessedDataDir() +
#'                    "/aberrant_splicing/datasets/savedObjects/raw-{dataset}/delta_theta.h5"`'
#'  type: script
#'--- 

saveRDS(snakemake, snakemake@log$snakemake)
suppressPackageStartupMessages(library(FRASER))

# input
dataset <- snakemake@wildcards$dataset
fds_in  <- file.path(dirname(snakemake@input$counting_done), "fds-object.RDS")
params  <- snakemake@config$aberrantSplicing
BPPARAM  <- MulticoreParam(snakemake@threads)

# Set number of threads including for DelayedArray operations
register(BPPARAM)
DelayedArray::setAutoBPPARAM(BPPARAM)

# Force writing HDF5 files
options(FRASER.maxSamplesNoHDF5=-1)
options(FRASER.maxJunctionsNoHDF5=-1)

# load FraserDataSet object
fds <- loadFraserDataSet(file=fds_in)

# Calculating PSI values
fds <- calculatePSIValues(fds)

# FRASER object after PSI value calculation
fds <- saveFraserDataSet(fds)
