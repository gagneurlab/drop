#'---
#' title: Filter and clean dataset
#' author: Christian Mertes
#' wb:
#'  log:
#'    - snakemake: '`sm str(tmp_dir / "AS" / "{dataset}" / "03_filter.Rds")`'
#'  params:
#'   - setup: '`sm cfg.AS.getWorkdir() + "/config.R"`'
#'   - workingDirIn: '`sm cfg.getProcessedDataDir() + "/aberrant_splicing/datasets/fromBam/"`'
#'   - workingDirOut: '`sm cfg.getProcessedDataDir() + "/aberrant_splicing/datasets/merged/"`'
#'   - exCountIDs: '`sm lambda w: sa.getIDsByGroup(w.dataset, assay="SPLICE_COUNT")`'
#'  input:
#'   - theta:  '`sm cfg.getProcessedDataDir()+
#'                  "/aberrant_splicing/datasets/fromBam/savedObjects/raw-{dataset}/theta.h5"`'
#'   - exCounts: '`sm lambda w: cfg.AS.getExternalCounts(w.dataset, "k_j_counts")`'
#'  output:
#'   - fds: '`sm cfg.getProcessedDataDir() +
#'                "/aberrant_splicing/datasets/merged/savedObjects/{dataset}/fds-object.RDS"`'
#'   - done: '`sm cfg.getProcessedDataDir() + 
#'                "/aberrant_splicing/datasets/merged/savedObjects/{dataset}/filter.done" `'
#'  threads: 3
#'  type: script
#'---

saveRDS(snakemake, snakemake@log$snakemake)
source(snakemake@params$setup, echo=FALSE)

opts_chunk$set(fig.width=12, fig.height=8)

# input
dataset    <- snakemake@wildcards$dataset
workingDirIn <- snakemake@params$workingDirIn
workingDirOut <- snakemake@params$workingDirOut
params     <- snakemake@config$aberrantSplicing
exCountIDs <- snakemake@params$exCountIDs
exCountFiles <- snakemake@input$exCounts
sample_anno_file <- snakemake@config$sampleAnnotation
minExpressionInOneSample <- params$minExpressionInOneSample
minDeltaPsi <- params$minDeltaPsi

fds <- loadFraserDataSet(dir=workingDirIn, name=paste0("raw-", dataset))

register(MulticoreParam(snakemake@threads))
# Limit number of threads for DelayedArray operations
setAutoBPPARAM(MulticoreParam(snakemake@threads))

# Add external data if provided by dataset
if(length(exCountIDs) > 0){
    message("create new merged fraser object")
    workingDir(fds) <- workingDirOut
    fds <- saveFraserDataSet(fds,dir = workingDirOut, name=paste0("raw-", dataset))

    for(resource in unique(exCountFiles)){
        exSampleIDs <- exCountIDs[exCountFiles == resource]
        exAnno <- fread(sample_anno_file, key="RNA_ID")[J(exSampleIDs)]
        setnames(exAnno, "RNA_ID", "sampleID")
        
        ctsNames <- c("k_j", "k_theta", "n_psi3", "n_psi5", "n_theta")
        ctsFiles <- paste0(dirname(resource), "/", ctsNames, "_counts.tsv.gz")
        
        fds <- mergeExternalData(fds=fds, countFiles=ctsFiles,
                sampleIDs=exSampleIDs, annotation=exAnno)
    }
} else {
    message("symLink fraser dir")
    file.symlink(paste0(workingDirIn, "savedObjects/","raw-", dataset),
                 paste0(workingDirOut, "savedObjects/","raw-", dataset))
    
    workingDir(fds) <- workingDirOut
    name(fds) <- paste0("raw-", dataset)
}

# filter for expression and write it out to disc.
# 
# TODO:   This will brake a rerun of step 01_5_countRNA_collect.R as it writes 
#         out the rawCountsJ and rawCountsSS file including the external samples. 
# 
fds <- filterExpressionAndVariability(fds, 
        minExpressionInOneSample = minExpressionInOneSample,
        minDeltaPsi = minDeltaPsi,
        filter=FALSE)

devNull <- saveFraserDataSet(fds,dir = workingDirOut)

# Keep junctions that pass filter
name(fds) <- dataset
if (params$filter == TRUE) {
    filtered <- mcols(fds, type="j")[,"passed"]
    fds <- fds[filtered,]
    message(paste("filtered to", nrow(fds), "junctions"))
}

seqlevels(fds) <- seqlevelsInUse(fds)
colData(fds)$sampleID <- as.character(colData(fds)$sampleID)
fds <- saveFraserDataSet(fds,dir = workingDirOut)
file.create(snakemake@output$done)
