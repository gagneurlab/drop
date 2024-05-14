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
#'   - splice_metrics: '`sm expand(cfg.getProcessedDataDir() +
#'                  "/aberrant_splicing/datasets/savedObjects/raw-local-{dataset}/{type}.h5", type=cfg.AS.getPsiTypeAssay(), allow_missing=True)`'
#'   - exCounts: '`sm lambda w: cfg.AS.getExternalCounts(w.dataset, "k_j_counts")`'
#'  output:
#'   - fds: '`sm cfg.getProcessedDataDir() +
#'                  "/aberrant_splicing/datasets/savedObjects/{dataset}/fds-object.RDS"`'
#'   - done: '`sm expand(cfg.getProcessedDataDir() +
#'                  "/aberrant_splicing/datasets/savedObjects/{dataset}/filter_{version}.done", version=cfg.AS.get("FRASER_version"), allow_missing=True)`'
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
quantile <- params$quantileForFiltering
quantileMinExpression <- params$quantileMinExpression
minDeltaPsi <- params$minDeltaPsi
filterOnJaccard <- (params$FRASER_version == "FRASER2")

fds <- loadFraserDataSet(dir=workingDir, name=paste0("raw-local-", dataset))

register(MulticoreParam(snakemake@threads))
# Limit number of threads for DelayedArray operations
setAutoBPPARAM(MulticoreParam(snakemake@threads))

# Add external data if provided by dataset
if(length(exCountIDs) > 0){
    message("create new merged fraser object")
    fds <- saveFraserDataSet(fds,dir = workingDir, name=paste0("raw-", dataset))

    for(resource in unique(exCountFiles)){
        exSampleIDs <- exCountIDs[exCountFiles == resource]
        exAnno <- fread(sample_anno_file, key="RNA_ID")[J(exSampleIDs)]
        exAnno[, strand:=exAnno$STRAND]
        exAnno$strand<- sapply(exAnno[, strand], switch, 'no' = 0L, 'unstranded' = 0L, 
                                   'yes' = 1L, 'stranded' = 1L, 'reverse' = 2L, -1L)
        setnames(exAnno, "RNA_ID", "sampleID")
        
        ctsNames <- c("k_j", "k_theta", "n_psi3", "n_psi5", "n_theta")
        ctsFiles <- paste0(dirname(resource), "/", ctsNames, "_counts.tsv.gz")
        
        # Merging external counts restricts the junctions to those that 
        # are only present in both the counted (fromBam) junctions AND the 
        # junctions from the external counts.
        fds <- mergeExternalData(fds=fds, countFiles=ctsFiles,
                sampleIDs=exSampleIDs, annotation=exAnno)
        fds@colData$isExternal <- as.factor(!is.na(fds@colData$SPLICE_COUNTS_DIR))
    }
} else {
    message("symLink fraser dir")
    file.symlink(paste0(workingDir, "savedObjects/","raw-local-", dataset),
                 paste0(workingDir, "savedObjects/","raw-", dataset))
    
    fds@colData$isExternal <- as.factor(FALSE)
    workingDir(fds) <- workingDir
    name(fds) <- paste0("raw-", dataset)
}

# filter for expression and write it out to disc.
fds <- filterExpressionAndVariability(fds, 
        minExpressionInOneSample = minExpressionInOneSample,
        quantile=quantile,
        quantileMinExpression=quantileMinExpression,
        minDeltaPsi = minDeltaPsi,
        filterOnJaccard=filterOnJaccard,
        filter=FALSE)

devNull <- saveFraserDataSet(fds,dir = workingDir)

# Keep junctions that pass filter
name(fds) <- dataset
if (params$filter == TRUE) {
    filtered <- mcols(fds, type="j")[,"passed"]
    fds <- fds[filtered,]
    message(paste("filtered to", nrow(fds), "junctions"))
}

seqlevels(fds) <- seqlevelsInUse(fds)
colData(fds)$sampleID <- as.character(colData(fds)$sampleID)
fds <- saveFraserDataSet(fds,dir = workingDir)

# remove previous filter.done files and create new one
outdir <- dirname(snakemake@output$done)
prevFilterFiles <- grep("filter(.*)done", list.files(outdir), value=TRUE)
unlink(file.path(outdir, prevFilterFiles))
file.create(snakemake@output$done)
