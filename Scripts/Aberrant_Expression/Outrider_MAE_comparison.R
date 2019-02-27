#'---
#' title: OUTRIDER and MAE results comparison
#' author: vyepez
#' wb:
#'  input:
#'   - ods_ss: '`sm expand(config["PROC_RESULTS"] + "/{annotation}/outrider/ss/ods.Rds", annotation=config["ANNOTATIONS"])`'
#'   - ods_ns: '`sm expand(config["PROC_RESULTS"] + "/{annotation}/outrider/ns/ods.Rds", annotation=config["ANNOTATIONS"])`'
#'   - mae_res: '`sm config["PROC_RESULTS"] + "/mae/MAE_results.Rds"`'
#'  output:
#' output: 
#'   html_document:
#'    code_folding: show
#'    code_download: TRUE
#'---

#+ echo=F
saveRDS(snakemake, 'tmp/mae_out.Rds')
# snakemake <- readRDS('tmp/mae_out.Rds')
suppressPackageStartupMessages({
    library(data.table)
    library(magrittr)
    library(tidyr)
    library(ggplot2)
    library(OUTRIDER)
})

#'
ods_ss <- readRDS(snakemake@input$ods_ss)
dt1 <- data.table(aberrant(ods_ss) %>% colSums())
dt1[, RNA_ID := colnames(ods_ss)]

ods_ns <- readRDS(snakemake@input$ods_ns)
dt2 <- data.table(aberrant(ods_ns) %>% colSums())
dt2[, RNA_ID := colnames(ods_ns)]

dt_ab <- rbind(dt1, dt2)
setnames(dt_ab, "V1", "outrider_aberrant")

mae_res <- readRDS(snakemake@input$mae_res)

mae_res <- separate(mae_res, "sample", into = c("EXOME_ID", "RNA_ID"), sep = "-")
mae_res_sample <- mae_res[, .(mae_events = .N), by = RNA_ID]
mae_samples <- unique(mae_res[,.(EXOME_ID, RNA_ID)])
mae_samples[duplicated(mae_samples$EXOME_ID)]

ab_dt <- merge(dt_ab, mae_res_sample)
ggplot(ab_dt, aes(outrider_aberrant, mae_events)) + geom_point() + theme_bw()

sa <- fread("../sample_annotation/Data/sample_annotation.tsv")
sa[RNA_ID %in% ab_dt[mae_events > 500, RNA_ID]]


