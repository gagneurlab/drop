#'---
#' title: Calculate P values
#' author: Christian Mertes
#' wb:
#'  params:
#'   - tmpdir: '`sm drop.getMethodPath(METHOD, "tmp_dir")`'
#'   - workingDir: '`sm parser.getProcDataDir() + "/aberrant_splicing/datasets/"`'
#'  threads: 20
#'  input:
#'   - fdsin:  '`sm parser.getProcDataDir() + 
#'                  "/aberrant_splicing/datasets/savedObjects/{dataset}/" +
#'                  "predictedMeans_psiSite.h5"`'
#'  output:
#'   - fdsout: '`sm parser.getProcDataDir() + 
#'                  "/aberrant_splicing/datasets/savedObjects/{dataset}/" +
#'                  "padjBetaBinomial_psiSite.h5"`'
#'  type: script
#'---

saveRDS(snakemake, file.path(snakemake@params$tmpdir, "FRASER_06.snakemake"))
# snakemake <- readRDS(".drop/tmp/AS/FRASER_06.snakemake")

source("Scripts/_helpers/config.R")

dataset    <- snakemake@wildcards$dataset
fdsFile    <- snakemake@input$fdsin
workingDir <- snakemake@params$workingDir

register(MulticoreParam(snakemake@threads))
# Limit number of threads for DelayedArray operations
setAutoBPPARAM(MulticoreParam(snakemake@threads))

# Load Zscores data
fds <- loadFraserDataSet(dir=workingDir, name=dataset)

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

