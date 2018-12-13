#'---
#' title: Count reads
#' author: Michaela Muller
#' wb:
#'  input:
#'   - sample_bam: '`sm config["RAW_DATA"] + "/{sampleID}/RNAout/paired-endout/stdFilenames/{sampleID}.bam"`'
#'   - features: '`sm "resources/exons_op.Rds"`'
#'  output:
#'   - counts: '`sm config["PROC_DATA"] + "/{sampleID}_counts_unspec.RDS"`'
#'   - wBhtml: 'Output/html/{sampleID}-counts.html'
#'  type: noindex
#'---

saveRDS(snakemake, "tmp/counts.snakemake")

message(paste("input:", snakemake@input$features))
message(paste("output:", snakemake@output$counts))

bam_file <- Rsamtools::BamFile(snakemake@input$sample_bam, yieldSize = 2e6)
feature_regions <- readRDS(snakemake@input$features)

# Subset for Testing
param <- Rsamtools::ScanBamParam(which = GenomicRanges::GRanges("chr21", IRanges::IRanges(1, 536870912)))
feature_regions <- feature_regions[GenomicRanges::seqnames(feature_regions) == "chr21"
                                   & GenomicRanges::end(feature_regions) <= 12000000,]

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
    , inter.feature = T   	# TRUE, reads mapping to multiple features are dropped
    , param = param
)
saveRDS(se, snakemake@output$counts)
message("done")

print(format(Sys.time()- starttime)) # time taken
print(sum(SummarizedExperiment::assay(se))) # total counts


