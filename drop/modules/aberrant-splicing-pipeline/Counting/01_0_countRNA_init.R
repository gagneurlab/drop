#'---
#' title: Initialize Counting
#' author: Luise Schuller
#' wb:
#'  log:
#'    - snakemake: '`sm str(tmp_dir / "AS" / "{dataset}" / "01_0_init.Rds")`'
#'  params:
#'   - setup: '`sm cfg.AS.getWorkdir() + "/config.R"`'
#'   - workingDir: '`sm cfg.getProcessedDataDir() + "/aberrant_splicing/datasets/"`'
#'  input:
#'    - colData: '`sm cfg.getProcessedDataDir() + 
#'                    "/aberrant_splicing/annotations/{dataset}.tsv"`'
#'  output:
#'   - fdsobj:  '`sm cfg.getProcessedDataDir() + 
#'                   "/aberrant_splicing/datasets/savedObjects/raw-local-{dataset}/fds-object.RDS"`'
#'   - done_fds: '`sm cfg.getProcessedDataDir() + 
#'                "/aberrant_splicing/datasets/cache/raw-local-{dataset}/fds.done" `'
#'  type: script
#'---

saveRDS(snakemake, snakemake@log$snakemake)
source(snakemake@params$setup, echo=FALSE)

dataset    <- snakemake@wildcards$dataset
colDataFile <- snakemake@input$colData
workingDir <- snakemake@params$workingDir
params <- snakemake@config$aberrantSplicing

# Create initial FRASER object
col_data <- fread(colDataFile)

col_data$strand <- 0L

fds <- FraserDataSet(colData = col_data,
                     workingDir = workingDir,
                     name       = paste0("raw-local-", dataset))

# Add paired end and strand specificity to the fds
pairedEnd(fds) <- colData(fds)$PAIRED_END
strandSpecific(fds) <- colData(fds)$STRAND

# Save initial FRASER dataset
fds <- saveFraserDataSet(fds)

message(date(), ": FRASER object initialized for ", dataset)

file.create(snakemake@output$done_fds)
