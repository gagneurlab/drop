#'---
#' title: Count reads
#' author: Michaela Mueller
#' wb:
#'  log:
#'    snakemake: '`sm str(tmp_dir / "AE" / "{annotation}" / "counts" / "{sampleID}.Rds")`'
#'  params:
#'   - COUNT_PARAMS: '`sm lambda w: cfg.AE.getCountParams(w.sampleID)`'
#'  input:
#'   - sample_bam: '`sm lambda w: sa.getFilePath(w.sampleID, file_type="RNA_BAM_FILE") `'
#'   - count_ranges: '`sm cfg.getProcessedDataDir() + "/aberrant_expression/{annotation}/count_ranges.Rds" `'
#'   - input_params: '`sm cfg.getProcessedDataDir() + "/aberrant_expression/{annotation}/params/counts/{sampleID}_countParams.csv" `'
#'  output:
#'   - counts: '`sm cfg.getProcessedDataDir() + "/aberrant_expression/{annotation}/counts/{sampleID,[^/]+}.Rds"`'
#'  type: script
#'  threads: 1
#'---

saveRDS(snakemake, snakemake@log$snakemake)

suppressPackageStartupMessages({
  library(data.table)
  library(Rsamtools)
  library(BiocParallel)
  library(GenomicAlignments)
})

# Get strand specific information from sample annotation
sampleID <- snakemake@wildcards$sampleID

count_params <- snakemake@params$COUNT_PARAMS
strand <- tolower(count_params$STRAND)
count_mode <- count_params$COUNT_MODE
paired_end <- as.logical(count_params$PAIRED_END)
overlap <- as.logical(count_params$COUNT_OVERLAPS)
inter_feature <- ! overlap # inter_feature = FALSE does not allow overlaps

# infer preprocessing and strand info
preprocess_reads <- NULL
if (strand == "yes") {
  strand_spec <- T
} else if (strand == "no"){
  strand_spec <- F
} else if (strand == "reverse") {
  # set preprocess function for later
  preprocess_reads <- invertStrand
  strand_spec <- T
} else {
  stop(paste("invalid strand information", strand))
}

# read files
bam_file <- BamFile(snakemake@input$sample_bam, yieldSize = snakemake@config$aberrantExpression$yieldSize)
count_ranges <- readRDS(snakemake@input$count_ranges)
# set chromosome style
seqlevelsStyle(count_ranges) <- seqlevelsStyle(bam_file)

# show info
message(paste("input:", snakemake@input$sample_bam))
message(paste("output:", snakemake@output$counts))
message(paste('\tcount mode:', count_mode, sep = "\t"))
message(paste('\tpaired end:', paired_end, sep = "\t"))
message(paste('\tinter.feature:', inter_feature, sep = "\t"))
message(paste('\tstrand:', strand, sep = "\t"))
message(paste('\tstrand specific:', strand_spec, sep = "\t"))
message(paste(seqlevels(count_ranges), collapse = ' '))

# start counting
message("\ncounting")
starttime <- Sys.time()
se <- summarizeOverlaps(
  count_ranges
    , bam_file
    , mode = count_mode
    , singleEnd = !paired_end
    , ignore.strand = !strand_spec  # FALSE if done strand specifically
    , fragments = F
    , count.mapped.reads = T
    , inter.feature = inter_feature # TRUE: reads mapping to multiple features are dropped
    , preprocess.reads = preprocess_reads
    , BPPARAM = MulticoreParam(snakemake@threads)
)
colnames(se) <- sampleID
saveRDS(se, snakemake@output$counts)
message("done")

print(format(Sys.time()- starttime)) # print time taken
print(sum(assay(se))) # total counts

