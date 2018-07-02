#'---
#' title: Count reads of samples on extracted features
#' author: Michaela Müller
#' wb:
#'  input: 
#'  - sample_bam: '`sm config["RAW_DATA"] + "/{sampleID}R/RNAout/paired-endout/stdFilenames/{sampleID}R.bam", sampleID=config["SAMPLE_IDS"]`'
#'  output:
#'  - counts: '`sm config["DATA_DIR"] + "/total_counts.RDS"`'
#'---

saveRDS(snakemake, "tmp/count_all.RDS")

cores <- min(BiocParallel::bpworkers(), 30)
BiocParallel::register(BiocParallel::MulticoreParam(cores))

single_counts <- BiocParallel::bplapply(snakemake@input$counts, readRDS)
total_counts <- do.call(SummarizedExperiment::cbind, single_counts)

saveRDS(total_counts, snakemake@output$merged_counts)

