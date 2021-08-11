#'---
#' title: Calculate P values
#' author: Christian Mertes
#' wb:
#'  log:
#'    - snakemake: '`sm str(tmp_dir / "AS" / "{dataset}" / "06_stats.Rds")`'
#'  threads: 20
#'  input:
#'   - fdsin:  '`sm cfg.getProcessedDataDir() + 
#'                  "/aberrant_splicing/datasets/savedObjects/{dataset}/" +
#'                  "predictedMeans_theta.h5"`'
#'  output:
#'   - fdsout: '`sm cfg.getProcessedDataDir() + 
#'                  "/aberrant_splicing/datasets/savedObjects/{dataset}/" +
#'                  "padjBetaBinomial_theta.h5"`'
#'  type: script
#'---

saveRDS(snakemake, snakemake@log$snakemake)
suppressPackageStartupMessages(library(FRASER))

#+ input
dataset  <- snakemake@wildcards$dataset
fdsFile  <- snakemake@input$fdsin
BPPARAM  <- MulticoreParam(snakemake@threads)

# Set number of threads including for DelayedArray operations
register(BPPARAM)
DelayedArray::setAutoBPPARAM(BPPARAM)

# Force writing HDF5 files
options(FRASER.maxSamplesNoHDF5=-1)
options(FRASER.maxJunctionsNoHDF5=-1)

# Load Zscores data
fds <- loadFraserDataSet(file=fdsFile)

# Calculate stats
for (type in psiTypes) {
    # Zscores
    fds <- calculateZscore(fds, type=type)
    # Pvalues
    fds <- calculatePvalues(fds, type=type)
    # Adjust Pvalues
    fds <- calculatePadjValues(fds, type=type)
}

fds <- saveFraserDataSet(fds)

