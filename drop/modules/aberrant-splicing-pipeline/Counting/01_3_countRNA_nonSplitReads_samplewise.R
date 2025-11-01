#'---
#' title: Nonsplit Counts
#' author: Luise Schuller
#' wb:
#'  log:
#'    - snakemake: '`sm str(tmp_dir / "AS" / "{dataset}" / "nonsplitReads" / "{sample_id}.Rds")`'
#'  params:
#'   - setup: '`sm cfg.AS.getWorkdir() + "/config.R"`'
#'   - workingDir: '`sm cfg.getProcessedDataDir() + "/aberrant_splicing/datasets"`'
#'  input:
#'   - spliceSites: '`sm cfg.getProcessedDataDir() + 
#'                   "/aberrant_splicing/datasets/cache/raw-local-{dataset}/spliceSites_splitCounts.rds"`'
#'  output:
#'   - done_sample_nonSplitCounts : '`sm cfg.getProcessedDataDir() + 
#'                   "/aberrant_splicing/datasets/cache/raw-local-{dataset}/sample_tmp/nonSplitCounts/sample_{sample_id}.done"`' 
#'  threads: 3
#'  resources:
#'    - mem_mb: 3000
#'  type: script
#'---

saveRDS(snakemake, snakemake@log$snakemake)
source(snakemake@params$setup, echo=FALSE)

dataset    <- snakemake@wildcards$dataset
colDataFile <- snakemake@input$colData
workingDir <- snakemake@params$workingDir
params <- snakemake@config$aberrantSplicing

# Read FRASER object
fds <- loadFraserDataSet(dir=workingDir, name=paste0("raw-local-", dataset))

# Get sample id from wildcard
sample_id <- snakemake@wildcards[["sample_id"]]


# Read splice site coordinates from RDS
spliceSiteCoords <- readRDS(snakemake@input$spliceSites)

# Copy to local tmp
bamFile <- bamFile(fds[, samples(fds) == sample_id])[[1]]
tmp_file <- file.path(snakemake@resources$tmpdir, paste0('nonsplit_counts_', sample_id, '.bam'))

message('Copy ', bamFile, ' to ', tmp_file)
file.copy(bamFile, tmp_file)
file.copy(paste0(bamFile, '.bai'), paste0(tmp_file, '.bai'))
bamFile(fds[, samples(fds) == sample_id]) <- tmp_file

# Count nonSplitReads for given sample id
sample_result <- countNonSplicedReads(sample_id,
                                      splitCountRanges = NULL,
                                      fds = fds,
                                      NcpuPerSample = snakemake@threads,
                                      minAnchor=5,
                                      recount=params$recount,
                                      spliceSiteCoords=spliceSiteCoords,
                                      longRead=params$longRead)

message(date(), ": ", dataset, ", ", sample_id,
        " no. splice junctions (non split counts) = ", length(sample_result))

# Reset bam file
bamFile(fds[, samples(fds) == sample_id]) <- bamFile
file.remove(tmp_file)

file.create(snakemake@output$done_sample_nonSplitCounts)
