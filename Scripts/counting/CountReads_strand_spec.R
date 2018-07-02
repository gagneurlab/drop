#'---
#' title: Count reads of samples on extracted features
#' author: Michaela Müller
#' wb:
#'  input: 
#'  - sample_bam: '`sm config["RAW_DATA"] + "/{sampleID}R/RNAout/paired-endout/stdFilenames/{sampleID}R.bam"`'
#'  - features: '`sm "resources/exons_op.Rds"`'
#'  output:
#'  - counts: '`sm config["PROC_DATA"] + "/{sampleID}_counts.RDS"`'
#'  - wBhtml: 'Output/html/{sampleID}/counts.html'
#'  type: noindex
#'---

saveRDS(snakemake, "tmp/counts.snakemake")

message(paste("input:", snakemake@input$features))
message(paste("output:", snakemake@output$counts))

bam_file <- Rsamtools::BamFile(snakemake@input$sample_bam,yieldSize=2e6)
feature_regions <- readRDS(snakemake@input$features)

# Subset for Testing
param <- Rsamtools::ScanBamParam(which = GenomicRanges::GRanges("chr21", IRanges::IRanges(1, 536870912)))
feature_regions <- feature_regions[GenomicRanges::seqnames(feature_regions) == "chr21"
                                   & GenomicRanges::end(feature_regions) <= 100,]

message("counting")
starttime <- Sys.time()
se <- GenomicAlignments::summarizeOverlaps(
    feature_regions
    , bam_file
    , mode = "Union"
    , singleEnd = F
    , ignore.strand = F
    , fragments = F
    , count.mapped.reads = T
    , inter.feature = F
    # , param = param
)
saveRDS(se, snakemake@output$counts)
message("done")

print(format(Sys.time()- starttime)) # time taken
print(sum(SummarizedExperiment::assay(se))) # total counts


