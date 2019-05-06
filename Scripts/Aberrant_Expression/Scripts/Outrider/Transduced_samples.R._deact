#'---
#' title: Transduced gene RNAs
#' author: vyepez
#' wb:
#'  input:
#'   - ods_all: '`sm config["PROC_RESULTS"] + "/v29_overlap/outrider/all/ods.Rds"`'
#'   - res_all: '`sm config["PROC_RESULTS"] + "/v29_overlap/outrider/all/OUTRIDER_results.tsv"`'
#'  output:
#' output:
#'   html_document:
#'    code_folding: show
#'    code_download: TRUE
#'---

#+ echo=F
saveRDS(snakemake, "tmp/transduced_outrider.snakemake")
# snakemake <- readRDS("tmp/transduced_outrider.snakemake")

suppressPackageStartupMessages({
    library(OUTRIDER)
    library(ggplot2)
    library(data.table)
    library(dplyr)
    library(cowplot)
})

source("Scripts/_functions/gene_annotation/add_gene_info_cols.R")

#' ## Results
ods <- readRDS(snakemake@input$ods_all)
ods <- readRDS("/s/project/genetic_diagnosis/processed_results/v29_overlap/outrider/all/ods.Rds")
res <- fread(snakemake@input$res_all)
res <- fread("/s/project/genetic_diagnosis/processed_results/v29_overlap/outrider/all/OUTRIDER_results.tsv")
sa <- fread("../sample_annotation/Data/sample_annotation.tsv")

#' Annotated transduced samples
unique(sa[!is.na(TRANSDUCED_GENE),.(FIBROBLAST_ID, RNA_ID, TRANSDUCED_GENE)])

res <- left_join(res, sa[,.(RNA_ID, TRANSDUCED_GENE)], by = c("sampleID" = "RNA_ID")) %>% as.data.table()
res <- res[!is.na(TRANSDUCED_GENE)]
res <- cbind(res[,.(geneID, sampleID, TRANSDUCED_GENE)], res[, -c("geneID", "sampleID", "TRANSDUCED_GENE")])
res[, tp_trans := any(TRANSDUCED_GENE %in% geneID), by = sampleID]
if(res[sampleID == 'MUC1389' & geneID == 'C1QBP', .N] > 0) res[sampleID=='MUC1389', tp_trans:=T]
if(res[sampleID == 'MUC1430' & geneID == 'NSUN3', .N] > 0) res[sampleID=='MUC1430', tp_trans:=T]

res[, FIBROBLAST_ID := paste(FIBROBLAST_ID, TRANSDUCED_GENE, sep = "-")]

#' Transduced samples
unique(res[, .(FIBROBLAST_ID, sampleID, TRANSDUCED_GENE)])
sort(table(unique(res[, .(sampleID, TRANSDUCED_GENE)])[,TRANSDUCED_GENE]))


plot(sort(unique(res[, .(sampleID, AberrantBySample)])[,AberrantBySample]), xlab = 'Rank', ylab = 'Number of outlier genes')
unique(res[AberrantBySample>20, .(sampleID, FIBROBLAST_ID)])


unique(res[tp_trans == T, .(sampleID, FIBROBLAST_ID)])

unique(res[tp_trans == F, .(sampleID, FIBROBLAST_ID)])

OUTRIDER::results(ods["SSBP1", "99590R-T-SSBP1"], all = TRUE)
plotExpressionRank(ods, "SSBP1")

OUTRIDER::results(ods["C1QBP", "MUC1389"], all = TRUE)
plotExpressionRank(ods, "C1QBP")

plotExpressionRank(ods, "TXN2")
plotExpressionRank(ods, "FLAD1")

# TODO: Volcano plots of transduced samples
