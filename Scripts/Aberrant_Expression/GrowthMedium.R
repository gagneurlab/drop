#'---
#' title: Compare Expression across Growth Media
#' author: mumichae
#' wb:
#'  input:
#'   - res_signif_all: '`sm config["PROC_RESULTS"] + "/mae/MAE_results.Rds"`'
#'   - res_signif_rare: '`sm config["PROC_RESULTS"] + "/mae/MAE_results_rare.Rds"`'
#' output: 
#'   html_document:
#'    code_folding: show
#'    code_download: TRUE
#'---

#+ echo=F
saveRDS(snakemake, 'tmp/medium.Rds')
# snakemake <- readRDS('tmp/medium.Rds')
suppressPackageStartupMessages({
    library(data.table)
    library(magrittr)
    library(ggplot2)
    library(scales)
    library(cowplot)
    library(ggpubr)
    library(plotly)
    library(tidyr)
})

