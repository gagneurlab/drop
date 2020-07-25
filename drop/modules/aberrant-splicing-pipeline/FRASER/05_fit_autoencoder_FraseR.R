#'---
#' title: Fitting the autoencoder
#' author: Christian Mertes
#' wb:
#'  log:
#'    - snakemake: '`sm str(tmp_dir / "AS" / "{dataset}" / "05_fit.Rds")`'
#'  params:
#'   - setup: '`sm cfg.AS.getWorkdir() + "/config.R"`'
#'   - workingDir: '`sm cfg.getProcessedDataDir() + "/aberrant_splicing/datasets/"`'
#'  threads: 20
#'  input:
#'   - hyper: '`sm cfg.getProcessedDataDir() + 
#'                "/aberrant_splicing/datasets/savedObjects/{dataset}/hyper.done" `'
#'  output:
#'   - fdsout: '`sm cfg.getProcessedDataDir() + 
#'                  "/aberrant_splicing/datasets/savedObjects/{dataset}/predictedMeans_psiSite.h5"`'
#'  type: script
#'---

saveRDS(snakemake, snakemake@log$snakemake)
source(snakemake@params$setup, echo=FALSE)

dataset    <- snakemake@wildcards$dataset
workingDir <- snakemake@params$workingDir

register(MulticoreParam(snakemake@threads))
# Limit number of threads for DelayedArray operations
setAutoBPPARAM(MulticoreParam(snakemake@threads))

fds <- loadFraserDataSet(dir=workingDir, name=dataset)

# Fit autoencoder
# run it for every type
implementation <- snakemake@config$aberrantSplicing$implementation

for(type in psiTypes){
    currentType(fds) <- type
    q <- bestQ(fds, type)
    fds <- fit(fds, q=q, type=type, verbose=TRUE, iterations=15, 
               implementation=implementation)
    fds <- saveFraserDataSet(fds)
}

