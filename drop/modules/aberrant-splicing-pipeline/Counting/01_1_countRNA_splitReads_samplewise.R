#'---
#' title: Count Split Reads
#' author: Luise Schuller
#' wb:
#'  log:
#'    - snakemake: '`sm str(tmp_dir / "AS" / "splitReads" / "{sampleID}.Rds")`'
#'  params:
#'    - COUNT_PARAMS: '`sm lambda w: cfg.AE.getCountParams(w.sampleID)`'
#'  input:
#'    - sample_bam: '`sm lambda w: sa.getFilePath(w.sampleID, file_type="RNA_BAM_FILE") `'
#'  output:
#'    - sample_splitCounts: '`sm cfg.getProcessedDataDir() + 
#'                   "/aberrant_splicing/datasets/cache/splitCounts-{sampleID}.RDS"`'
#'  threads: 3
#'  type: script
#'---

saveRDS(snakemake, snakemake@log$snakemake)
suppressPackageStartupMessages({
  library(FRASER)
  library(BSgenome)
})

sampleID   <- snakemake@wildcards$sampleID
sample_bam <- snakemake@input$sample_bam
fds_init   <- snakemake@input$fds_init
out_rds    <- snakemake@output$sample_splitCounts
params <- snakemake@config$aberrantSplicing
genomeAssembly <- snakemake@config$genomeAssembly
strand <- snakemake@params$COUNT_PARAMS$STRAND
paired <- snakemake@params$COUNT_PARAMS$PAIRED_END
coldata <- fread(snakemake@config$sampleAnnotation)[RNA_ID == sampleID]


# Get strand specificity of sample ('' does not exists in switch -> x_...)
strand <- switch(paste0("x_", tolower(strand)),
        'x_' = 0L,
        'x_no' = 0L,
        'x_unstranded' = 0L,
        'x_yes' = 1L,
        'x_stranded' = 1L,
        'x_reverse' = 2L,
        stop("Strand parameter not known! It was: '", strand, "'."))

# If data is not strand specific, add genome info
genome <- NULL
if(strand == 0){
  genome <- getBSgenome(genomeAssembly)
}

# Count splitReads for a given sample id
sample_result <- countSplitReads(
    sampleID = sampleID,
    NcpuPerSample = snakemake@threads,
    genome = genome,
    recount = TRUE,
    keepNonStandardChromosomes=params$keepNonStandardChrs,
    bamfile = sample_bam,
    pairedend = paired,
    strandmode = strand,
    cacheFile = out_rds,
    scanbamparam = ScanBamParam(mapqFilter=0), 
    coldata = coldata)

message(date(), ": ", sampleID," no. splice junctions (split counts) = ", 
    length(sample_result))
