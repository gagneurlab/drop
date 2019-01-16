#'---
#' title: Count reads
#' author: Michaela Mueller
#' wb:
#'  input:
#'   - sample_bam: '`sm config["RAW_DATA"] + "/{sampleID}/RNAout/paired-endout/stdFilenames/{sampleID}.bam"`'
#'   - features: '`sm config["PROC_RESULTS"] + "/{annotation}/counts/exons_by_gene_op.Rds"`'
#'  output:
#'   - counts: '`sm config["PROC_RESULTS"] + "/counts/overlap/{annotation}/{sampleID}_counts.Rds"`'
#'  type: script
#'---

#source(".wBuild/wBuildParser.R")
#parseWBHeader("Scripts/counting/countReads.R")
saveRDS(snakemake, "tmp/counts.snakemake")

message(paste("input:", snakemake@input$features))
message(paste("output:", snakemake@output$counts))

bam_file <- Rsamtools::BamFile(snakemake@input$sample_bam, yieldSize = 2e6)
feature_regions <- readRDS(snakemake@input$features)

message("counting")
starttime <- Sys.time()
se <- GenomicAlignments::summarizeOverlaps(
    feature_regions
    , bam_file
    , mode = 'IntersectionStrict'
    , singleEnd = F
    , ignore.strand = F  # FALSE if done strand specifically
    , fragments = F
    , count.mapped.reads = T
    , inter.feature = F   	# TRUE, reads mapping to multiple features are dropped
)
saveRDS(se, snakemake@output$counts)
message("done")

print(format(Sys.time()- starttime)) # time taken
print(sum(SummarizedExperiment::assay(se))) # total counts


