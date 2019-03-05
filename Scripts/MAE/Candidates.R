#'---
#' title: MAE Candidates
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
saveRDS(snakemake, 'tmp/candidates.Rds')
# snakemake <- readRDS('tmp/candidates.Rds')
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

res <- readRDS(snakemake@input$res_signif_all)
res_rare <- readRDS(snakemake@input$res_signif_rare)
DT::datatable(res_rare[!is.na(HANS_CLASS) & is.na(KNOWN_MUTATION)], caption = "Unsolved mitodisease MAE results", style = 'bootstrap')

total_res_MUC1404 <- readRDS("/s/project/genetic_diagnosis/processed_results/mae/samples/65990-MUC1404_res.Rds")
total_res <- total_res[MAX_AF < 0.2]
total_res_MUC1404[, fc := 2^(log2FC)]
total_res_MUC1404[, label := 'unsignificant']
total_res_MUC1404[padj < 0.5 & alt_freq > .8, label := 'monoallelic']
total_res_MUC1404[label == 'monoallelic' & MAX_AF <= 10e-3, label := 'rare monoallelic']
g <- ggplot(total_res_MUC1404, aes(alt_cov + ref_cov, fc, col = label)) +
    geom_point(aes(text = hgncid)) +
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
    labs(title = 'RNA: MUC1404, Exome: 65990, Patient: #80256', y = 'fold change of (ALT+1)/(REF+1)') +
    scale_color_manual(values = c('cornflowerblue', 'firebrick', 'grey80')) +
    theme(legend.position = 'bottom')
g
ggplotly(g)

