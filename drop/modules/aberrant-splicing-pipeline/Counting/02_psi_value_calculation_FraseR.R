#'---
#' title: Calculate PSI values
#' author: Christian Mertes
#' wb:
#'  log:
#'   - snakemake: '`sm str(tmp_dir / "AS" / "{dataset}" / "02_PSIcalc.Rds")`'
#'  params:
#'   - setup: '`sm cfg.AS.getWorkdir() + "/config.R"`'
#'   - workingDir: '`sm cfg.getProcessedDataDir() + "/aberrant_splicing/datasets/"`'
#'  threads: 30
#'  input:
#'   - counting_done: '`sm cfg.getProcessedDataDir() + 
#'                "/aberrant_splicing/datasets/savedObjects/raw-local-{dataset}/counting.done" `'
#'  output:
#'  - splice_metrics: '`sm expand(cfg.getProcessedDataDir() +
#'                    "/aberrant_splicing/datasets/savedObjects/raw-local-{dataset}/{type}.h5", type=cfg.AS.getPsiTypeAssay(), allow_missing=True)`'
#'  type: script
#'--- 

saveRDS(snakemake, snakemake@log$snakemake)
source(snakemake@params$setup, echo=FALSE)

dataset    <- snakemake@wildcards$dataset
workingDir <- snakemake@params$workingDir
params <- snakemake@config$aberrantSplicing

register(MulticoreParam(snakemake@threads))
# Limit number of threads for DelayedArray operations
setAutoBPPARAM(MulticoreParam(snakemake@threads))

fds <- loadFraserDataSet(dir=workingDir, name=paste0("raw-local-", dataset))

# Calculating PSI values
fds <- calculatePSIValues(fds)

# FRASER object after PSI value calculation
fds <- saveFraserDataSet(fds)
