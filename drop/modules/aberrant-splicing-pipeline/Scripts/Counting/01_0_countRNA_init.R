#'---
#' title: Initialize Counting
#' author: Luise Schuller
#' wb:
#'  params:
#'   - tmpdir: '`sm drop.getMethodPath(METHOD, "tmp_dir")`'
#'   - workingDir: '`sm parser.getProcDataDir() + "/aberrant_splicing/datasets"`'
#'  input:
#'    - colData: '`sm parser.getProcDataDir() + 
#'                    "/aberrant_splicing/annotations/{dataset}.tsv"`'
#'  output:
#'   - fdsobj:  '`sm parser.getProcDataDir() + 
#'                   "/aberrant_splicing/datasets/savedObjects/raw-{dataset}/fds-object.RDS"`'
#'   - done_fds: '`sm parser.getProcDataDir() + 
#'                "/aberrant_splicing/datasets/cache/raw-{dataset}/fds.done" `'
#'  type: script
#'---
saveRDS(snakemake, file.path(snakemake@params$tmpdir, "FRASER_01_0.snakemake"))
# snakemake <- readRDS(".drop/tmp/AS/FRASER_01_0.snakemake")

source("Scripts/_helpers/config.R")

dataset    <- snakemake@wildcards$dataset
colDataFile <- snakemake@input$colData
workingDir <- snakemake@params$workingDir
params <- snakemake@config$aberrantSplicing

# Create initial FRASER object
col_data <- fread(colDataFile)

fds <- FraserDataSet(colData = col_data,
                     workingDir = workingDir,
                     name       = paste0("raw-", dataset))

# Add paired end and strand specificity to the fds
pairedEnd(fds) <- colData(fds)$PAIRED_END
strandSpecific(fds) <- 'no'
if(uniqueN(colData(fds)$STRAND) == 1){
  strandSpecific(fds) <- unique(colData(fds)$STRAND)
} 

# Save initial FRASER dataset
fds <- saveFraserDataSet(fds)

message(date(), ": FRASER object initialized for ", dataset)

file.create(snakemake@output$done_fds)
