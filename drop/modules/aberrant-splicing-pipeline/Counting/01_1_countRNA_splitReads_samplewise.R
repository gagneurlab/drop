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

# Count splitReads for a given sample id
sample_result <- countSplitReads(sampleID = sample_id, 
                                 fds = fds,
                                 NcpuPerSample = snakemake@threads,
                                 recount = params$recount,
                                 keepNonStandardChromosomes = params$keepNonStandardChrs,
                                 genome = genome)

message(date(), ": ", dataset, ", ", sample_id,
        " no. splice junctions (split counts) = ", length(sample_result))

file.create(snakemake@output$done_sample_splitCounts)
