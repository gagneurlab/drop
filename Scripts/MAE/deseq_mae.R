#'---
#' title: Get MAE results
#' author: vyepez
#' wb:
#'  input:
#'   - mae_counts: '`sm config["PROC_DATA"] + "/mae/{vcf}-{rna}.Rds"`'
#'   - vcf_uniqs: '`sm config["RAW_DATA"] + "/{vcf}/exomicout/paired-endout/processedData/vep_anno_{vcf}_uniq_dt.Rds"`'
#'  output:
#'   - mae_res: '`sm config["PROC_RESULTS"] + "/mae/samples/{vcf}-{rna}_res.Rds"`'
#'  type: script
#'---

saveRDS(snakemake, 'tmp/res_mae.Rds')
# snakemake <- readRDS(snakemake, 'tmp/res_mae.Rds')

suppressPackageStartupMessages({
    devtools::load_all("../mae/")
    library(dplyr)
    library(data.table)
    library(magrittr)
    library(tidyr)
})

mae_raw <- readRDS(snakemake@input$mae_counts)

# Function from MAE pkg
rmae <- run_deseq_all_mae(mae_raw)

# Make an aux column for merging
rmae[, aux := paste(chr, pos, REF, ALT, sep = "-")]

# Read vcf annotated file with 1 row / variant
vt <- readRDS(snakemake@input$vcf_uniqs)
vt[, chr := paste0("chr", chr)]
vt[, aux := paste(chr, pos, ref, alt, sep = "-")]

# Merge results
rt <- left_join(rmae, vt[,.(hgncid, mstype, noccds, sift1, pph1, gnomAD_AF, MAX_AF, AF, gnomAD_NFE_AF, gnomAD_AFR_AF, gnomAD_EAS_AF, rsid, pubmed, aux)], by = "aux") %>% as.data.table
rt[, aux := NULL]

# Save results
saveRDS(rt, snakemake@output$mae_res)
