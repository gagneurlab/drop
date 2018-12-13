#'---
#' title: Count reads of samples on extracted features
#' author: Michaela Müller
#' wb:
#'  input: 
#'  - counts: '`sm expand(config["PROC_DATA"] + "/{sampleID}_counts_unspec.RDS", sampleID=config["SAMPLE_IDS"])`'
#'  output:
#'  - tot_counts: '`sm config["PROC_DATA"] + "/total_counts.RDS"`'
#'---

saveRDS(snakemake, "tmp/count_all.RDS")

cores <- min(BiocParallel::bpworkers(), 30)
BiocParallel::register(BiocParallel::MulticoreParam(cores))

list = paste0(PROC_DATA, samples, "_counts_unspec.RDS")
single_counts <- BiocParallel::bplapply(list, readRDS)
single_counts <- BiocParallel::bplapply(snakemake@input$counts, readRDS)
total_counts <- do.call(SummarizedExperiment::cbind, single_counts)

saveRDS(total_counts, snakemake@output$merged_counts)


samples <- c("102860R", "102864R", "102874R")
    