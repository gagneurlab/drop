#'---
#' title: Initialize Counting
#' author: Luise Schuller
#' wb:
#'  log:
#'    - snakemake: '`sm str(tmp_dir / "AS" / "{dataset}" / "01_0_init.Rds")`'
#'  input:
#'    - colData: '`sm cfg.getProcessedDataDir() + 
#'                    "/aberrant_splicing/annotations/{dataset}.tsv"`'
#'  output:
#'   - fds_init:  '`sm cfg.getProcessedDataDir() + 
#'                   "/aberrant_splicing/datasets/savedObjects/raw-{dataset}/fds-init.done"`'
#'  type: script
#'---

saveRDS(snakemake, snakemake@log$snakemake)
suppressPackageStartupMessages(library(FRASER))

# input
dataset     <- snakemake@wildcards$dataset
colDataFile <- snakemake@input$colData
fds_init    <- snakemake@output$fds_init
params      <- snakemake@config$aberrantSplicing

# Create initial FRASER object
col_data <- fread(colDataFile)

fds <- FraserDataSet(colData = col_data,
                     workingDir = dirname(dirname(dirname(fds_init))),
                     name       = paste0("raw-", dataset))

# Add strand specificity to the fds
strandSpecific(fds) <- 'no'
if(uniqueN(colData(fds)$STRAND) == 1){
  strandSpecific(fds) <- unique(colData(fds)$STRAND)
} 

# Save initial FRASER dataset
fds <- saveFraserDataSet(fds)

message(date(), ": FRASER object initialized for ", dataset)

file.create(fds_init)
