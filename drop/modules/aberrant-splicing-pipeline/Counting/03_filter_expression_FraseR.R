#'---
#' title: Filter and clean dataset
#' author: Christian Mertes
#' wb:
#'  log:
#'    - snakemake: '`sm str(tmp_dir / "AS" / "{dataset}" / "03_filter.Rds")`'
#'  params:
#'   - workingDir: '`sm cfg.getProcessedDataDir() + "/aberrant_splicing/datasets/"`'
#'  input:
#'   - setup: '`sm cfg.AS.getWorkdir() + "/config.R"`'
#'   - theta:  '`sm cfg.getProcessedDataDir()+
#'                  "/aberrant_splicing/datasets/savedObjects/raw-{dataset}/theta.h5"`'
#'  output:
#'   - fds: '`sm cfg.getProcessedDataDir() +
#'                "/aberrant_splicing/datasets/savedObjects/{dataset}/fds-object.RDS"`'
#'   - done: '`sm cfg.getProcessedDataDir() + 
#'                "/aberrant_splicing/datasets/savedObjects/{dataset}/filter.done" `'
#'  threads: 3
#'  type: script
#'---

saveRDS(snakemake, snakemake@log$snakemake)
source(snakemake@input$setup, echo=FALSE)
#source(snakemake@input$spliceTypeSetup, echo=FALSE)
#source(snakemake@input$addAnnotation)

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

# Add the junction annotations to the fds
#message("03: load db for annotation")
#txdb <- loadDb(snakemake@input$txdb)
#seqlevelsStyle(txdb) <- seqlevelsStyle(fds)
#fds <- createFDSAnnotations(fds, txdb)
#message("03: save object after annotation")

fds <- saveFraserDataSet(fds)
file.create(snakemake@output$done)
