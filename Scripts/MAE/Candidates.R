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

plot_volcano <- function(mae_res_dt, sample = NULL, interactive = T) {
    mae_res_dt <- mae_res_dt[MAX_AF < 0.2]
    mae_res_dt[, fc := 2^(log2FC)]
    mae_res_dt[, label := 'unsignificant']
    mae_res_dt[padj < 0.5 & alt_freq > .8, label := 'monoallelic']
    mae_res_dt[label == 'monoallelic' & MAX_AF <= 10e-3, label := 'rare monoallelic']
    
    g <- ggplot(mae_res_dt, aes(alt_cov + ref_cov, fc, col = label)) +
        geom_point(aes(text = hgncid)) +
        labs(title = 'RNA: MUC1404, Exome: 65990, Patient: #80256', y = 'fold change of (ALT+1)/(REF+1)') +
        scale_color_manual(values = c('dodgerblue', 'orangered', 'grey80')) +
        theme(legend.position = 'bottom')
    
    if (interactive) {
        ggplotly(g + scale_x_log10() + scale_y_log10())
    }
    else {
        g + scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
            scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))
    }
}

plot_volcano(total_res_MUC1404, interactive = F)

