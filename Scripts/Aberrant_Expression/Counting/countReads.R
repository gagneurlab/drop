#'---
#' title: Count reads
#' author: Michaela Mueller
#' wb:
#'  input:
#'   - sample_bam: '`sm config["RAW_DATA"] + "/{sampleID}/RNAout/paired-endout/stdFilenames/{sampleID}.bam"`'
#'   - features: '`sm config["PROC_RESULTS"] + "/{annotation}/counts/exons_by_gene_op.Rds"`'
#'  output:
#'   - counts: '`sm config["PROC_RESULTS"] + "/{annotation}/counts/{sampleID}_counts.Rds"`'
#'  type: script
#'---

#source(".wBuild/wBuildParser.R")
#parseWBHeader("Scripts/counting/countReads.R")
saveRDS(snakemake, "tmp/counts.snakemake")
# snakemake <- readRDS("tmp/counts.snakemake")
suppressPackageStartupMessages({
    library(data.table)
})

# import count settings from config
anno <- snakemake@wildcards$annotation
dt <- data.table(annotation = unlist(snakemake@config["ANNOTATIONS"]),
                 inter_feature = as.logical(unlist(snakemake@config["INTER_FEATURE"])))

message(paste("input:", snakemake@input$features))
message(paste("output:", snakemake@output$counts))
message("settings:")
message(dt)
message(paste('anno:', anno))
message(paste('inter.feature:', dt[annotation == anno, inter_feature]))


# read files
bam_file <- Rsamtools::BamFile(snakemake@input$sample_bam, yieldSize = 2e6)
feature_regions <- readRDS(snakemake@input$features)

# start counting
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
    , inter.feature = dt[annotation == anno, inter_feature] # TRUE, reads mapping to multiple features are dropped
)
saveRDS(se, snakemake@output$counts)
message("done")

print(format(Sys.time()- starttime)) # time taken
print(sum(SummarizedExperiment::assay(se))) # total counts


