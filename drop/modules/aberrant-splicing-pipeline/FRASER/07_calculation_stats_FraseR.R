#'---
#' title: Calculate P values
#' author: Christian Mertes
#' wb:
#'  log:
#'    - snakemake: '`sm str(tmp_dir / "AS" / "{dataset}--{annotation}" / "07_stats.Rds")`'
#'  params:
#'   - setup: '`sm cfg.AS.getWorkdir() + "/config.R"`'
#'   - workingDir: '`sm cfg.getProcessedResultsDir() + "/aberrant_splicing/datasets/"`'
#'   - ids: '`sm lambda w: sa.getIDsByGroup(w.dataset, assay="RNA")`'
#'   - parse_subsets_for_FDR: '`sm str(projectDir / ".drop" / "helpers" / "parse_subsets_for_FDR.R")`'
#'   - genes_to_test: '`sm cfg.AS.get("genesToTest")`'
#'  threads: 20
#'  input:
#'   - fdsin:  '`sm expand(cfg.getProcessedResultsDir() +
#'                 "/aberrant_splicing/datasets/savedObjects/{dataset}--{annotation}/" +
#'                  "predictedMeans_{type}.h5", type=cfg.AS.getPsiTypeAssay(), allow_missing=True)`'
#'  output:
#'   - fdsout:  '`sm expand(cfg.getProcessedResultsDir() +
#'                 "/aberrant_splicing/datasets/savedObjects/{dataset}--{annotation}/" +
#'                  "padjBetaBinomial_{type}.h5", type=cfg.AS.getPsiTypeAssay(), allow_missing=True)`'
#'  type: script
#'---

saveRDS(snakemake, snakemake@log$snakemake)
source(snakemake@params$setup, echo=FALSE)

annotation <- snakemake@wildcards$annotation
dataset    <- snakemake@wildcards$dataset
fdsFile    <- snakemake@input$fdsin
workingDir <- snakemake@params$workingDir

register(MulticoreParam(snakemake@threads))
# Limit number of threads for DelayedArray operations
setAutoBPPARAM(MulticoreParam(snakemake@threads))

# read in subsets from sample anno if present (returns NULL if not present)
source(snakemake@params$parse_subsets_for_FDR)
fraser_sample_ids <- snakemake@params$ids
subsets <- parse_subsets_for_FDR(snakemake@params$genes_to_test,
                                 sampleIDs=fraser_sample_ids)

# Load FRASER data
fds <- loadFraserDataSet(dir=workingDir, name=paste(dataset, annotation, sep = '--'))

# Calculate stats
for (type in psiTypes) {
    # Pvalues
    fds <- calculatePvalues(fds, type=type)
    # Adjust Pvalues
    fds <- calculatePadjValues(fds, type=type, subsets=subsets)
}

fds <- saveFraserDataSet(fds)

# remove .h5 files from previous runs with other FRASER version
fdsDir <- dirname(snakemake@output$fdsout[1])
pvalFiles <- grep("p(.*)BetaBinomial_(.*).h5", 
                  list.files(fdsDir), 
                  value=TRUE)
for(type in psiTypesNotUsed){
    pvalFilesType <- grep(type, pvalFiles, value=TRUE)
    for(pFile in pvalFilesType){
        unlink(file.path(fdsDir, pFile))
    }
}
