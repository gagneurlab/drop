#'---
#' title: Run MAE for a sample
#' author: mumichae
#' wb:
#'  input:
#'   - vcf: '`sm config["RAW_DATA"] + "/{vcf}/exomicout/paired-endout/stdFilenames/{vcf}.bam"`'
#'   - rna: '`sm config["RAW_DATA"] + "/{rna}/RNAout/paired-endout/stdFilenames/{rna}.bam"`'
#'  output:
#'   - mae: '`sm config["PROC_DATA"] + "/mae/{vcf}-{rna}.Rds"`'
#'  threads: 5
#'  type: script
#'---

saveRDS(snakemake, 'tmp/mae.Rds')
suppressPackageStartupMessages({
    devtools::load_all("../mae/")
})

vcfs <- snakemake@input$vcf
rnas <- snakemake@input$rna

BPPARAM = MulticoreParam(snakemake@threads, snakemake@threads, progressbar=TRUE)
gr <- countMAEReads(vcfs[1], rnas[1], BPPARAM=BPPARAM)

saveRDS(gr, snakemake@output$mae)
