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
#'  threads: 6
#'  resources:
#'    - mem_mb: '`sm lambda wildcards, threads, attempt: threads * 2000 + 10000 * attempt`'
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
yield_size <- snakemake@config$aberrantExpression$yieldSize
threads <- snakemake@threads

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
bam_file <- BamFile(snakemake@input$sample_bam, yieldSize = yield_size)
count_ranges <- readRDS(snakemake@input$count_ranges)
# set chromosome style
seqlevelsStyle(count_ranges) <- seqlevelsStyle(bam_file)

# run it in parallel across all chromosomes
gene_seqnames = as.character(sapply(seqnames(count_ranges), runValue))

# show info
message(paste("input:", snakemake@input$sample_bam))
message(paste("output:", snakemake@output$counts))
message(paste('\tcount mode:', count_mode, sep = "\t"))
message(paste('\tpaired end:', paired_end, sep = "\t"))
message(paste('\tinter.feature:', inter_feature, sep = "\t"))
message(paste('\tstrand:', strand, sep = "\t"))
message(paste('\tstrand specific:', strand_spec, sep = "\t"))
message(paste(seqlevels(count_ranges), collapse = ' '))

# function to count per chromosome
count_per_chromosome <- function(i, gene_seqnames, count_ranges, bam_file, yield_size,
                                 count_mode, paired_end, strand_spec, preprocess_reads) {
  genes_to_count <- gene_seqnames == i
  message(paste("counting seqname '", i, "' with '", sum(genes_to_count), "' genes.", sep = ""))
  
  # subset the count_ranges for the current seqname
  current_ranges <- count_ranges[genes_to_count]
  
  # create a new BamFile for the current seqname
  current_bam_file <- BamFile(path(bam_file), yieldSize=yield_size)

  # summarize overlaps
  summarizeOverlaps(
    current_ranges,
    current_bam_file,
    mode = count_mode,
    singleEnd = !paired_end,
    ignore.strand = !strand_spec,
    fragments = FALSE,
    count.mapped.reads = TRUE,
    preprocess.reads = preprocess_reads,
    inter.feature = inter_feature,
    param = ScanBamParam(which = GRanges(i, IRanges(1, seqlengths(bam_file)[i])))
  )
}

# run counting in parallel across all chromosomes
message("\nstarting counting expression ...")
starttime <- Sys.time()

iterate <- seqlevels(count_ranges)
bpparam <- MulticoreParam(workers = threads, tasks = length(iterate))
counts <- bplapply(iterate, count_per_chromosome, 
          gene_seqnames, count_ranges, bam_file, yield_size,
          count_mode, paired_end, strand_spec, preprocess_reads,
          BPPARAM = bpparam
)

# merge SE objects - concatenate by rows (genes) across chromosomes
se <- do.call(rbind, counts)

colnames(se) <- sampleID
saveRDS(se, snakemake@output$counts)
message("done")

# print time and stats taken
print(paste("Time taken:", format(Sys.time()- starttime), "with", threads, "threads."))
print(paste("Total counts:", sum(assay(se))))

