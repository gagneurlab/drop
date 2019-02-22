#'---
#' title: Run MAE for a sample
#' author: mumichae
#' wb:
#'  input:
#'   - vcf: '`sm config["RAW_DATA"] + "/{vcf}/exomicout/paired-endout/stdFilenames/{vcf}.vcf.gz"`'
#'  output:
#'   - vcf: '`sm config["RAW_DATA"] + "/{vcf}/exomicout/paired-endout/processedData/vep_anno_{vcf}.vcf.gz"`'
#'  type: script
#'---

saveRDS(snakemake, 'tmp/variant_annotation.Rds')
# snakemake <-  readRDS('tmp/variant_annotation.Rds')

suppressPackageStartupMessages({
    devtools::load_all("../mae/")
})

source("Scripts/_functions/annotation_with_vep.R")

vep_param <- get_vep_params(version=90, num_forks=10, vcfFile=snakemake@output$vcf)
resCall <- ensemblVEP(snakemake@input$vcf, vep_param)  # The vep_param already contains the output file
message(resCall)

# input_vcfs <- snakemake@input$vcf
# annot_vcfs <- snakemake@output$vcf
# 
# register(MulticoreParam(snakemake@threads))
# bplapply(seq_along(input_vcfs), function(i) {
#     vep_param <- get_vep_params(version=90, num_forks=10, vcfFile=annot_vcfs[[i]])
#     resCall <- ensemblVEP(input_vcfs[[i]], vep_param)  # The vep_param already contains the output file
#     message(paste(vcf_file, resCall, sep = '\n'))
# })
