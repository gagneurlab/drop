#'---
#' title: Nonsplit Counts
#' author: Luise Schuller
#' wb:
#'  log:
#'    - snakemake: '`sm str(tmp_dir / "AS" / "{dataset}" / "nonsplitReads" / "{sample_id}.Rds")`'
#'  input:
#'   - spliceSites: '`sm cfg.getProcessedDataDir() + 
#'                   "/aberrant_splicing/datasets/cache/raw-{dataset}/spliceSites_splitCounts.rds"`'
#'   - fds_init:    '`sm cfg.getProcessedDataDir() + 
#'                   "/aberrant_splicing/datasets/savedObjects/raw-{dataset}/fds-init.done"`'
#'  output:
#'   - nonSplitCounts_sample : '`sm cfg.getProcessedDataDir() + 
#'                   "/aberrant_splicing/datasets/cache/nonSplicedCounts/raw-{dataset}/nonSplicedCounts-{sample_id}.h5"`'
#'  threads: 3
#'  type: script
#'---

saveRDS(snakemake, snakemake@log$snakemake)
suppressPackageStartupMessages(library(FRASER))

dataset    <- snakemake@wildcards$dataset
sample_id  <- snakemake@wildcards$sample_id
fds_init   <- snakemake@input$fds_init
params     <- snakemake@config$aberrantSplicing

# Read FRASER object
fds <- loadFraserDataSet(file=file.path(dirname(fds_init), "fds-object.RDS"))

# Read splice site coordinates from RDS
spliceSiteCoords <- readRDS(snakemake@input$spliceSites)

# Count nonSplitReads for given sample id
sample_result <- countNonSplicedReads(sample_id,
        splitCountRanges = NULL,
        fds = fds,
        NcpuPerSample = snakemake@threads,
        minAnchor=5,
        # TODO should be TRUE in the future as snakemake is checking this
        recount=params$recount,
        spliceSiteCoords=spliceSiteCoords,
        longRead=params$longRead)

message(date(), ": ", dataset, ", ", sample_id,
        " no. splice junctions (non split counts) = ", length(sample_result))
