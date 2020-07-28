#'---
#' title: Collect all counts from FRASER Object
#' author: Michaela Mueller, vyepez
#' wb:
#'  params:
#'   - tmpdir: '`sm drop.getMethodPath(METHOD, "tmp_dir")`'
#'   - workingDir: '`sm parser.getProcDataDir() + "/aberrant_splicing/datasets"`'
#'  input:
#'   - gRangesSplitCounts: '`sm parser.getProcDataDir() + 
#'                          "/aberrant_splicing/datasets/cache/raw-{dataset}/gRanges_splitCounts.rds"`'
#'   - spliceSites: '`sm parser.getProcDataDir() + 
#'                   "/aberrant_splicing/datasets/cache/raw-{dataset}/spliceSites_splitCounts.rds"`'
#'   - counting_done: '`sm parser.getProcDataDir() + 
#'                "/aberrant_splicing/datasets/savedObjects/raw-{dataset}/counting.done" `'
#'  output:
#'    - split_counts: '`sm parser.getProcResultsDir() + "/exported_counts/{dataset}--{genomeAssembly}--{annotation}/"
#'                + "splitCounts.tsv.gz"`'
#'    - nonsplit_counts: '`sm parser.getProcResultsDir() + "/exported_counts/{dataset}--{genomeAssembly}--{annotation}/"
#'                + "spliceSiteOverlapCounts.tsv.gz"`'
#'  type: script
#'---
saveRDS(snakemake, file.path(snakemake@params$tmpdir, "export.snakemake"))
# snakemake <- readRDS(".drop/tmp/AS/export.snakemake")

source("Scripts/_helpers/config.R")

dataset    <- snakemake@wildcards$dataset
workingDir <- snakemake@params$workingDir


# Read FRASER object
fds <- loadFraserDataSet(dir=workingDir, name=paste0("raw-", dataset))
splitCounts_gRanges <- readRDS(snakemake@input$gRangesSplitCounts)
spliceSiteCoords <- readRDS(snakemake@input$spliceSites)

# obtain the split counts
splitCounts <- counts(fds, type="j")
gr_dt <- as.data.table(splitCounts_gRanges)
splitCounts <- cbind(gr_dt, as.matrix(splitCounts))
fwrite(splitCounts, file = snakemake@output$split_counts,
       quote = FALSE, row.names = FALSE, sep = '\t', compress = 'gzip')
  
# obtain the non split counts
nonSplitCounts <- counts(fds, type="ss")
grns_dt <- as.data.table(spliceSiteCoords)
nonSplitCounts <- cbind(grns_dt, as.matrix(nonSplitCounts))
fwrite(nonSplitCounts, file = snakemake@output$nonsplit_counts,
       quote = FALSE, row.names = FALSE, sep = '\t', compress = 'gzip')
