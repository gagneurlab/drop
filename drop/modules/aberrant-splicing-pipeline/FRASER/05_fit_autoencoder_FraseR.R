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
#'   - hyper: '`sm expand(cfg.getProcessedDataDir() + 
#'                "/aberrant_splicing/datasets/savedObjects/{dataset}/hyper_{version}.done", version=cfg.AS.get("FRASER_version"), allow_missing=True)`'
#'  output:
#'   - fdsout: '`sm expand(cfg.getProcessedDataDir() + 
#'                  "/aberrant_splicing/datasets/savedObjects/{dataset}/predictedMeans_{type}.h5", type=cfg.AS.getPsiTypeAssay(), allow_missing=True)`'
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
  verbose(fds) <- 3   # Add verbosity to the FRASER object
  fds <- fit(fds, q=q, type=type, iterations=15, implementation=implementation)
}

# remove .h5 files from previous runs with other FRASER version
fdsDir <- dirname(snakemake@output$fdsout[1])
for(type in psiTypesNotUsed){
  predMeansFile <- file.path(fdsDir, paste0("predictedMeans_", type, ".h5"))
  if(file.exists(predMeansFile)){
    unlink(predMeansFile)
    assay(fds, paste0("predictedMeans_", type), withDimnames=FALSE) <- NULL
  }
  deltaFiles <- file.path(fdsDir, paste0("delta_", type, ".h5"))
  if(file.exists(deltaFiles)){
    unlink(deltaFiles)
    assay(fds, paste0("delta_", type), withDimnames=FALSE) <- NULL
  }
  psiFiles <- file.path(fdsDir, paste0(type, ".h5"))
  if(file.exists(psiFiles)){
    unlink(psiFiles)
    assay(fds, type, withDimnames=FALSE) <- NULL
  }
  rawOtherFiles <- file.path(fdsDir, paste0("rawOtherCounts_", type, ".h5"))
  if(file.exists(rawOtherFiles)){
    unlink(rawOtherFiles)
    assay(fds, paste0("rawOtherCounts_", type), withDimnames=FALSE) <- NULL
  }
}

fds <- saveFraserDataSet(fds)
