#'---
#' title: Data table of annotated variants (1 per row)
#' author: mumichae
#' wb:
#'  input:
#'   - vcf: "{rawdata}/processedData/vep_anno_{vcf}.vcf.gz"
#'   - functions: "Scripts/_functions/annotation_with_vep.R"
#'  output:
#'   - vcf_dt: "{rawdata}/processedData/vep_anno_{vcf}_uniq_dt.Rds"
#'  type: script
#'---

saveRDS(snakemake, 'tmp/variant_annotation_dt.Rds')
# snakemake <-  readRDS('tmp/variant_annotation_dt.Rds')
source(snakemake@input$functions)

vcf_obj <- readVcf(snakemake@input$vcf, "hg19")
vcf_obj <- vcf_obj[rowRanges(vcf_obj)$QUAL >= 20]
vep_obj <- parseCSQToGRanges(vcf_obj, VCFRowID = rownames(vcf_obj))
vcf_dt <- combine_vcf_vep(sample = snakemake@wildcards$vcf, vcf_obj = vcf_obj, vep_obj = vep_obj)

saveRDS(vcf_dt, snakemake@output$vcf_dt)
