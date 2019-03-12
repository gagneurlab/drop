#'---
#' title: Run MAE for a sample
#' author: mumichae
#' wb:
#'  input:
#'   - vcf: '`sm standardFileNames("Data/helmholtz/{vcf}/exomicout/paired-endout/stdFilenames/{vcf}.vcf.gz")`'
#'   - rna: '`sm standardFileNames("Data/helmholtz/{rna}/RNAout/paired-endout/stdFilenames/{rna}.bam")`'
#'  output:
#'   - mae: '`sm config["PROC_DATA"] + "/mae/{vcf}-{rna}.Rds"`'
#'  threads: 1
#'  type: script
#'---

saveRDS(snakemake, 'tmp/mae.Rds')
# snakemake <-  readRDS('tmp/mae.Rds')

suppressPackageStartupMessages({
    devtools::load_all("../mae/")
    library(VariantAnnotation)
})

source("Scripts/_functions/filter_sets.R")
vcfs <- snakemake@input$vcf
rnas <- snakemake@input$rna

BPPARAM = MulticoreParam(snakemake@threads, snakemake@threads, progressbar=TRUE)
gr <- countMAEReads(vcfs, rnas, filter_function = filter_vcf_quality, BPPARAM=BPPARAM)  # already filters for quality

saveRDS(gr, snakemake@output$mae)
