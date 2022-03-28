#'---
#' title: Filter and clean dataset
#' author: Christian Mertes
#' wb:
#'  log:
#'    - snakemake: '`sm str(tmp_dir / "AS" / "{dataset}" / "03_filter.Rds")`'
#'  params:
#'   - setup: '`sm cfg.AS.getWorkdir() + "/config.R"`'
#'   - workingDir: '`sm cfg.getProcessedDataDir() + "/aberrant_splicing/datasets/"`'
#'   - exCountIDs: '`sm lambda w: sa.getIDsByGroup(w.dataset, assay="SPLICE_COUNT")`'
#'  input:
#'   - theta:  '`sm cfg.getProcessedDataDir()+
#'                  "/aberrant_splicing/datasets/savedObjects/raw-{dataset}/theta.h5"`'
#'   - exCounts: '`sm lambda w: cfg.AS.getExternalCounts(w.dataset, "k_j_counts")`'
#'  output:
#'   - fds: '`sm cfg.getProcessedDataDir() +
#'                "/aberrant_splicing/datasets/savedObjects/{dataset}/fds-object.RDS"`'
#'   - done: '`sm cfg.getProcessedDataDir() + 
#'                "/aberrant_splicing/datasets/savedObjects/{dataset}/filter.done" `'
#'  threads: 3
#'  type: script
#'---

saveRDS(snakemake, snakemake@log$snakemake)
source(snakemake@params$setup, echo=FALSE)

opts_chunk$set(fig.width=12, fig.height=8)

# input
dataset    <- snakemake@wildcards$dataset
workingDir <- snakemake@params$workingDir
params     <- snakemake@config$aberrantSplicing
exCountIDs <- snakemake@params$exCountIDs
exCountFiles <- snakemake@input$exCounts
sample_anno_file <- snakemake@config$sampleAnnotation
minExpressionInOneSample <- params$minExpressionInOneSample
minDeltaPsi <- params$minDeltaPsi

fds <- loadFraserDataSet(dir=workingDir, name=paste0("raw-", dataset))

register(MulticoreParam(snakemake@threads))
# Limit number of threads for DelayedArray operations
setAutoBPPARAM(MulticoreParam(snakemake@threads))

# Add external data if provided by dataset
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

# filter for expression and write it out to disc.
# 
# TODO:   This will brake a rerun of step 01_5_countRNA_collect.R as it writes 
#         out the rawCountsJ and rawCountsSS file including the external samples. 
# 
fds <- filterExpressionAndVariability(fds, 
        minExpressionInOneSample = minExpressionInOneSample,
        minDeltaPsi = minDeltaPsi,
        filter=FALSE)
fds <- saveFraserDataSet(fds)

# Keep junctions that pass filter
name(fds) <- dataset
if (params$filter == TRUE) {
    filtered <- mcols(fds, type="j")[,"passed"]
    fds <- fds[filtered,]
    message(paste("filtered to", nrow(fds), "junctions"))
}

seqlevels(fds) <- seqlevelsInUse(fds)
colData(fds)$sampleID <- as.character(colData(fds)$sampleID)
fds <- saveFraserDataSet(fds)
file.create(snakemake@output$done)
