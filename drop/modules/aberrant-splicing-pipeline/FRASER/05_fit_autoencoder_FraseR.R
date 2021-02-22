#'---
#' title: Fitting the autoencoder
#' author: Christian Mertes
#' wb:
#'  log:
#'    - snakemake: '`sm str(tmp_dir / "AS" / "{dataset}" / "05_fit.Rds")`'
#'  threads: 20
#'  input:
#'   - hyper: '`sm cfg.getProcessedDataDir() + 
#'                "/aberrant_splicing/datasets/savedObjects/{dataset}/hyper.done" `'
#'  output:
#'   - fdsout: '`sm cfg.getProcessedDataDir() + 
#'                  "/aberrant_splicing/datasets/savedObjects/{dataset}/predictedMeans_theta.h5"`'
#'  type: script
#'---

saveRDS(snakemake, snakemake@log$snakemake)
suppressPackageStartupMessages(library(FRASER))

# input
dataset <- snakemake@wildcards$dataset
fds_in  <- file.path(dirname(snakemake@input$hyper), "fds-object.RDS")
BPPARAM <- MulticoreParam(snakemake@threads)

# Set number of threads including for DelayedArray operations
register(BPPARAM)
DelayedArray::setAutoBPPARAM(BPPARAM)

# load FraserDataSet object
fds <- loadFraserDataSet(file=fds_in)

# Fit autoencoder
# run it for every type
implementation <- snakemake@config$aberrantSplicing$implementation

for(type in psiTypes){
    currentType(fds) <- type
    q <- bestQ(fds, type)
    verbose(fds) <- 3   # Add verbosity to the FRASER object
    fds <- fit(fds, q=q, type=type, iterations=15, implementation=implementation)
    fds <- saveFraserDataSet(fds)
}

