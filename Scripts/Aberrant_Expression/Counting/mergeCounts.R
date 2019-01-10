#'---
#' title: Merge the counts for all samples
#' author: Michaela Muller
#' wb:
#'  input: 
#'  - counts: '`sm expand(config["PROC_RESULTS"] + "/{{annotation}}/counts/{sampleID}_counts.Rds", sampleID=config["SAMPLE_IDS"])`'
#'  output:
#'  - merged_counts: '`sm config["PROC_RESULTS"] + "/{annotation}/counts/total_counts.Rds"`'
#'  threads: 60
#'  type: script
#'---

saveRDS(snakemake, "tmp/count_all.snakemake")
suppressPackageStartupMessages({
    library(BiocParallel)
    library(SummarizedExperiment)
})

register(MulticoreParam(snakemake@threads))

single_counts <- bplapply(snakemake@input$counts, readRDS)
message(paste("read", length(single_counts), "files"))

total_counts <- do.call(cbind, single_counts)
colnames(total_counts) <- gsub(".bam", "", colnames(total_counts))

saveRDS(total_counts, snakemake@output$merged_counts)
    