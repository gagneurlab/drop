#'---
#' title: OUTRIDER and MAE results comparison
#' author: vyepez
#' wb:
#'  input:
#'   - ods_ss: '`sm expand(config["PROC_RESULTS"] + "/{annotation}/outrider/fib_ss/ods.Rds", annotation=config["ANNOTATIONS"])`'
#'   - ods_ns: '`sm expand(config["PROC_RESULTS"] + "/{annotation}/outrider/fib_ns/ods.Rds", annotation=config["ANNOTATIONS"])`'
#'   - ods_res: '`sm expand(config["PROC_RESULTS"] + "/{annotation}/outrider/OUTRIDER_results.tsv", annotation=config["ANNOTATIONS"])`'
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

#' ## Read ods objects and results tables
# Get the number of significant outliers per sample
ods_ss <- readRDS(snakemake@input$ods_ss)
dt1 <- data.table(aberrant(ods_ss) %>% colSums())
dt1[, RNA_ID := colnames(ods_ss)]

ods_ns <- readRDS(snakemake@input$ods_ns)
dt2 <- data.table(aberrant(ods_ns) %>% colSums())
dt2[, RNA_ID := colnames(ods_ns)]

dt_ab <- rbind(dt1, dt2)
setnames(dt_ab, "V1", "outrider_aberrant")

mae_res <- readRDS(snakemake@input$mae_res)

mae_res_sample <- mae_res[, .(mae_events = .N), by = RNA_ID]
mae_samples <- unique(mae_res[,.(EXOME_ID, RNA_ID)])
mae_samples[duplicated(mae_samples$EXOME_ID)]

ab_dt <- merge(dt_ab, mae_res_sample)
ggplot(ab_dt, aes(outrider_aberrant, mae_events)) + geom_point() + theme_bw()

#' ## Mismatches
sa <- fread("../sample_annotation/Data/sample_annotation.tsv")
DT::datatable(sa[RNA_ID %in% ab_dt[mae_events > 500, RNA_ID]])


#' ## Venn Diagrams
ods_res <- fread(snakemake@input$ods_res)
ods_res[, gene_result := paste(sampleID, geneID, sep = '-')]
mae_res <- mae_res[!is.na(gene_name) & (MAX_AF <= 0.001 | is.na(MAX_AF))]
mae_res[, gene_result := paste(RNA_ID, gene_name, sep = '-')]
mae_res[, N_gr := .N, by = gene_result]
mae_res[, samples_in_common := RNA_ID %in% c(colnames(ods_ns), colnames(ods_ss))]
gplots::venn(list(MAE = unique(mae_res[samples_in_common == T, gene_result]), 
                  OUTRIDER = unique(ods_res$gene_result))
)
intersect(unique(mae_res[samples_in_common == T, gene_result]), unique(ods_res$gene_result))

