#'---
#' title: Create samples' table for MAE
#' author: vyepez
#' wb:
#'  input:
#'   - sample_anno: '/s/project/mitoMultiOmics/raw_data/sample_info/SAMPLE_ANNOTATION.tsv'
#'  output:
#'   - mae_sample_anno: '/s/project/genetic_diagnosis/processed_data/mae_files.tsv'
#' output:
#'   html_document:
#'    code_folding: show
#'    code_download: TRUE
#'---

#+ echo=F
saveRDS(snakemake, "tmp/sample_mae.snakemake")
# snakemake <- readRDS("tmp/sample_mae.snakemake")

suppressPackageStartupMessages({
    library(ggplot2)
    library(data.table)
})

sa <- fread(snakemake@input$sample_anno)

#'
mae_sa <- sa[,.(RNA_ID, EXOME_ID, FIBROBLAST_ID, BATCH, RNA_PERSON, TISSUE, COMMENT)]
dim(mae_sa)
mae_sa <- mae_sa[! grep("NHDF", FIBROBLAST_ID)]

#' From the following samples we have RNA, but no exome
DT::datatable(mae_sa[is.na(EXOME_ID)], style = 'bootstrap')

mae_sa[, RNA_file := paste0("/s/project/mitoMultiOmics/raw_data/helmholtz/", RNA_ID, "/RNAout/paired-endout/stdFilenames/", RNA_ID, ".bam")]
mae_sa[is.na(RNA_ID), RNA_file := NA]
mae_sa[, exome_vcf_file := paste0("/s/project/mitoMultiOmics/raw_data/helmholtz/", EXOME_ID, "/exomicout/paired-endout/stdFilenames/", EXOME_ID, ".vcf.gz")]
mae_sa[is.na(EXOME_ID), exome_vcf_file := NA]

mae_sa[, RNA_exists := file.exists(RNA_file)]
mae_sa[, vcf_exists := file.exists(exome_vcf_file)]

mae_sa[, both_files_exist := RNA_exists & vcf_exists]

#' From the following samples we have an RNA_ID, but no file
mae_sa[!is.na(RNA_ID) & RNA_exists == FALSE]

#' From the following samples we have an exome_ID, but no vcf file
DT::datatable(mae_sa[!is.na(EXOME_ID) & vcf_exists == FALSE])

ms <- mae_sa[both_files_exist == TRUE]

#' Number of samples with both RNA and VCF files
nrow(ms)

write.table(ms, snakemake@output$mae_sample_anno, quote = F, row.names = F, sep = "\t")