#'---
#' title: Collect all counts from FRASER Object
#' author: mumichae, vyepez
#' wb:
#'  log:
#'    - snakemake: '`sm str(tmp_dir / "AS" / "{dataset}" / "{genomeAssembly}--{annotation}_export.Rds")`'
#'  params:
#'   - setup: '`sm cfg.AS.getWorkdir() + "/config.R"`'
#'   - workingDir: '`sm cfg.getProcessedDataDir() + "/aberrant_splicing/datasets"`'
#'  input:
#'   - gRangesSplitCounts: '`sm cfg.getProcessedDataDir() + 
#'                          "/aberrant_splicing/datasets/cache/raw-{dataset}/gRanges_splitCounts.rds"`'
#'   - spliceSites: '`sm cfg.getProcessedDataDir() + 
#'                   "/aberrant_splicing/datasets/cache/raw-{dataset}/spliceSites_splitCounts.rds"`'
#'   - counting_done: '`sm cfg.getProcessedDataDir() + 
#'                "/aberrant_splicing/datasets/savedObjects/raw-{dataset}/counting.done" `'
#'  output:
#'    - split_counts: '`sm cfg.exportCounts.getFilePattern(str_=False) / "splitCounts.tsv.gz"`'
#'    - nonsplit_counts: '`sm cfg.exportCounts.getFilePattern(str_=False) / "spliceSiteOverlapCounts.tsv.gz"`'
#'  type: script
#'---

saveRDS(snakemake, snakemake@log$snakemake)
source(snakemake@params$setup, echo=FALSE)

dataset    <- snakemake@wildcards$dataset
workingDir <- snakemake@params$workingDir


# Read FRASER object
fds <- loadFraserDataSet(dir=workingDir, name=paste0("raw-", dataset))
splitCounts_gRanges <- readRDS(snakemake@input$gRangesSplitCounts)
spliceSiteCoords <- readRDS(snakemake@input$spliceSites)

# obtain the split counts
splitCounts <- counts(fds, type="j")
gr_dt <- as.data.table(splitCounts_gRanges)[, c(1:3,5)]
splitCounts <- cbind(gr_dt, as.matrix(splitCounts))
fwrite(splitCounts, file = snakemake@output$split_counts,
       quote = FALSE, row.names = FALSE, sep = '\t', compress = 'gzip')
  
# obtain the non split counts
nonSplitCounts <- counts(fds, type="ss")
grns_dt <- as.data.table(spliceSiteCoords)[, c(1:3,5)]
nonSplitCounts <- cbind(grns_dt, as.matrix(nonSplitCounts))
fwrite(nonSplitCounts, file = snakemake@output$nonsplit_counts,
       quote = FALSE, row.names = FALSE, sep = '\t', compress = 'gzip')
