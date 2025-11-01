#'---
#' title: Count Split Reads
#' author: Luise Schuller
#' wb:
#'  log:
#'    - snakemake: '`sm str(tmp_dir / "AS" / "{dataset}" / "splitReads" / "{sample_id}.Rds")`'
#'  params:
#'   - setup: '`sm cfg.AS.getWorkdir() + "/config.R"`'
#'   - workingDir: '`sm cfg.getProcessedDataDir() + "/aberrant_splicing/datasets/"`'
#'  input:
#'   - done_fds: '`sm cfg.getProcessedDataDir() + 
#'                "/aberrant_splicing/datasets/cache/raw-local-{dataset}/fds.done"`'
#'  output:
#'   - done_sample_splitCounts: '`sm cfg.getProcessedDataDir() + 
#'                "/aberrant_splicing/datasets/cache/raw-local-{dataset}"
#'                +"/sample_tmp/splitCounts/sample_{sample_id}.done"`'
#'  threads: 3
#'  resources:
#'    - mem_mb: 3000
#'  type: script
#'---

saveRDS(snakemake, snakemake@log$snakemake)
source(snakemake@params$setup, echo=FALSE)
library(BSgenome)

dataset    <- snakemake@wildcards$dataset
workingDir <- snakemake@params$workingDir
params <- snakemake@config$aberrantSplicing
genomeAssembly <- snakemake@config$genomeAssembly

# Read FRASER object
fds <- loadFraserDataSet(dir=workingDir, name=paste0("raw-local-", dataset))

# Get sample id from wildcard
sample_id <- snakemake@wildcards[["sample_id"]]

# If data is not strand specific, add genome info
genome <- NULL
if(strandSpecific(fds[, sample_id]) == 0){
  genome <- getBSgenome(genomeAssembly)
}

# Copy to local tmp
bamFile <- bamFile(fds[, samples(fds) == sample_id])[[1]]
tmp_file <- file.path(snakemake@resources$tmpdir, paste0('split_counts_', sample_id, '.bam'))

message('Copy ', bamFile, ' to ', tmp_file)
file.copy(bamFile, tmp_file)
file.copy(paste0(bamFile, '.bai'), paste0(tmp_file, '.bai'))
bamFile(fds[, samples(fds) == sample_id]) <- tmp_file

# Count splitReads for a given sample id
sample_result <- countSplitReads(sampleID = sample_id, 
                                 fds = fds,
                                 NcpuPerSample = snakemake@threads,
                                 recount = params$recount,
                                 keepNonStandardChromosomes = params$keepNonStandardChrs,
                                 genome = genome)

message(date(), ": ", dataset, ", ", sample_id,
        " no. splice junctions (split counts) = ", length(sample_result))

# Reset bam file
bamFile(fds[, samples(fds) == sample_id]) <- bamFile
file.remove(tmp_file)

file.create(snakemake@output$done_sample_splitCounts)
