source("Scripts/_functions/annotation_with_vep.R")

sample <- snakemake@wildcards$vcf

vcf_obj <- readVcf(vepVcfFile, "hg19")

# discard any mutations with less minQUAL
vcf_obj <- vcf_obj[rowRanges(vcf_obj)$QUAL >= minQUAL]
vep_obj <- parseCSQToGRanges(vcf_obj, VCFRowID = rownames(vcf_obj))
vcf_dt <- combine_vcf_vep(sample, vcf_obj, vep_obj)
