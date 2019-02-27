#'---
#' title: Annotate VCF file
#' author: mumichae
#' wb:
#'  input:
#'   - vcf: "{rawdata}/stdFilenames/{vcf}.vcf.gz"
#'  output:
#'   - vcf: "{rawdata}/processedData/vep_anno_{vcf}.vcf.gz"
#'   - vcf_html: "{rawdata}/processedData/vep_anno_{vcf}.vcf.gz_summary.html"
#'  type: script
#'  threads: 10
#'---

saveRDS(snakemake, 'tmp/variant_annotation.Rds')
# snakemake <-  readRDS('tmp/variant_annotation.Rds')

suppressPackageStartupMessages({
    devtools::load_all("../mae/")
})

source("Scripts/_functions/annotation_with_vep.R")

vep_param <- get_vep_params(version=94, num_forks=snakemake@threads, vcfFile=snakemake@output$vcf)
resCall <- ensemblVEP(snakemake@input$vcf, vep_param)  # The vep_param already contains the output files

if(resCall != 0){
    stop(resCall)
}
