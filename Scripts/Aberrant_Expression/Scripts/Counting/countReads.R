#'---
#' title: Count reads
#' author: Michaela Mueller
#' wb:
#'  input:
#'   - sample_bam: '`sm standardFileNames("Data/helmholtz/{sampleID}/RNAout/paired-endout/stdFilenames/{sampleID}.bam")`'
#'   - features: '`sm config["PROC_RESULTS"] + "/{annotation}/counts/exons_by_gene_op.Rds"`'
#'  output:
#'   - counts: '`sm config["PROC_RESULTS"] + "/{annotation}/counts/{sampleID,[^/]+}.Rds"`'
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
count_settings <- data.table(annotation = unlist(snakemake@config["ANNOTATIONS"]),
                             inter_feature = as.logical(unlist(snakemake@config["INTER_FEATURE"])))
inter_feature <- count_settings[annotation == anno, inter_feature]

# import sample annotation
sampleID <- snakemake@wildcards$sampleID
sample_anno <- fread(snakemake@config$SAMPLE_ANNOTATION)
# Get strand specific information from sample annotation
strand_spec <- sample_anno[RNA_ID == sampleID, COUNT_STRAND_SPECIFIC]

# show info
message(paste("input:", snakemake@input$features))
message(paste("output:", snakemake@output$counts))
message(paste('\tinter.feature:', inter_feature, sep = "\t"))
message(paste('\tstrand specific:', strand_spec, sep = "\t"))


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
    , ignore.strand = !strand_spec  # FALSE if done strand specifically
    , fragments = F
    , count.mapped.reads = T
    , inter.feature = inter_feature # TRUE, reads mapping to multiple features are dropped
)
saveRDS(se, snakemake@output$counts)
message("done")

print(format(Sys.time()- starttime)) # time taken
print(sum(SummarizedExperiment::assay(se))) # total counts


