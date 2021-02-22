#'---
#' title: Filter and clean dataset
#' author: Christian Mertes
#' wb:
#'  log:
#'    - snakemake: '`sm str(tmp_dir / "AS" / "{dataset}" / "03_filter.Rds")`'
#'  params:
#'   - exCountIDs: '`sm lambda w: sa.getIDsByGroup(w.dataset, assay="SPLICE_COUNT")`'
#'  input:
#'   - theta:  '`sm cfg.getProcessedDataDir()+
#'                  "/aberrant_splicing/datasets/savedObjects/raw-{dataset}/delta_theta.h5"`'
#'   - exCounts: '`sm lambda w: cfg.AS.getExternalCounts(w.dataset, "k_j_counts")`'
#'  output:
#'   - fds_k_j: '`sm cfg.getProcessedDataDir() +
#'                "/aberrant_splicing/datasets/savedObjects/{dataset}/rawCountsJ.h5"`'
#'   - fds_k_theta: '`sm cfg.getProcessedDataDir() +
#'                "/aberrant_splicing/datasets/savedObjects/{dataset}/rawCountsSS.h5"`'
#'  threads: 3
#'  type: script
#'---

saveRDS(snakemake, snakemake@log$snakemake)
suppressPackageStartupMessages(library(FRASER))
opts_chunk$set(fig.width=12, fig.height=8)

# input
dataset    <- snakemake@wildcards$dataset
fdsIn      <- snakemake@input$theta
params     <- snakemake@config$aberrantSplicing
exCountIDs <- snakemake@params$exCountIDs
exCountFiles <- snakemake@input$exCounts
sample_anno_file <- snakemake@config$sampleAnnotation
minExpressionInOneSample <- params$minExpressionInOneSample
minDeltaPsi <- params$minDeltaPsi
BPPARAM  <- MulticoreParam(snakemake@threads)

# Set number of threads including for DelayedArray operations
register(BPPARAM)
DelayedArray::setAutoBPPARAM(BPPARAM)

# Load FraserDataSet object
fds <- loadFraserDataSet(file=fdsIn)

# Apply filter
fds <- filterExpressionAndVariability(fds, 
        minExpressionInOneSample = minExpressionInOneSample,
        minDeltaPsi = minDeltaPsi,
        filter=FALSE)
devNull <- saveFraserDataSet(fds)

# TODO move it before applying filter to have it included in the summary
# Add external data if provided by dataset
dontWriteHDF5(fds) <- TRUE
if(length(exCountIDs) > 0){
    for(resource in unique(exCountFiles)){
        exSampleIDs <- exCountIDs[exCountFiles == resource]
        exAnno <- fread(sample_anno_file, key="RNA_ID")[J(exSampleIDs)]
        setnames(exAnno, "RNA_ID", "sampleID")
        
        ctsNames <- c("k_j", "k_theta", "n_psi3", "n_psi5", "n_theta")
        ctsFiles <- paste0(dirname(resource), "/", ctsNames, "_counts.tsv.gz")
        fds <- mergeExternalData(fds=fds, countFiles=ctsFiles,
                sampleIDs=exSampleIDs, annotation=exAnno)
    }
}

fds <- filterExpressionAndVariability(fds, 
        minExpressionInOneSample = minExpressionInOneSample,
        minDeltaPsi = minDeltaPsi,
        filter=params$filter)

# save filtered data into new dataset object
dontWriteHDF5(fds) <- FALSE
name(fds) <- dataset
fds <- saveFraserDataSet(fds, rewrite=TRUE)
