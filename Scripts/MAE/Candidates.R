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


get_total_mae_result <- function(sa, rna_id = NULL, exome_id = NULL, dir = "/s/project/genetic_diagnosis/processed_results/mae/samples") {
    if (!is.null(rna_id)){
        ids <- sa[RNA_ID == rna_id, .(RNA_ID, EXOME_ID)]
    } else if (!is.null(exome_id)) {
        ids <- sa[EXOME_ID == exome_id, .(RNA_ID, EXOME_ID)]
    } else {
        stop('either RNA or Exome ID have to be given')
    }
    if (nrow(ids) != 1)
        message(paste('number of RNA and exome IDs for', rna_id, 'and', exome_id, 'is not 1'))
    readRDS(file.path(dir, paste0(ids$EXOME_ID, '-', ids$RNA_ID, '_res.Rds')))
    # 65990-MUC1404_res.Rds
}

plot_volcano <- function(mae_res_dt, interactive = T, patient_id = NULL) {
    
    # mae_res_dt <- mae_res_dt[(MAX_AF < 0.2 | is.na(MAX_AF))]
    mae_res_dt[, fc := 2^(log2FC)]
    mae_res_dt[, label := 'unsignificant']
    mae_res_dt[padj < 0.5 & alt_freq > .8, label := 'monoallelic']
    mae_res_dt[label == 'monoallelic' & (MAX_AF <= 10e-3 | is.na(MAX_AF)), label := 'rare monoallelic']
    
    rna_id <- strsplit(unique(mae_res_dt$sample), '-')[[1]][1]
    exome_id <- strsplit(unique(mae_res_dt$sample), '-')[[1]][2]
    title <- paste0('RNA: ', rna_id, ', Exome: ', exome_id)
    if (!is.null(patient_id))
        title <- paste0(title, ', Patient: ', patient_id)
    
    g <- ggplot(mae_res_dt, aes(alt_cov + ref_cov, fc, col = label)) +
        geom_point(aes(text = hgncid)) +
        labs(title = title, y = 'fold change of (ALT+1)/(REF+1)') +
        scale_color_manual(values = c('dodgerblue', 'orangered', 'grey80')) +
        theme(legend.position = 'bottom')
    
    if (interactive) {
        ggplotly(g + scale_x_log10() + scale_y_log10())
    }
    else {
        g + scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
            scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
            grids()
    }
}


sa <- fread(snakemake@config$SAMPLE_ANNOTATION)
#+ volcanos, fig.height=8, fig.width=10
total_res_MUC1404 <- readRDS("/s/project/genetic_diagnosis/processed_results/mae/samples/65990-MUC1404_res.Rds")
plot_volcano(total_res_MUC1404, interactive = T, patient_id = '#80256')
plot_volcano(get_total_mae_result(sa, rna_id = '103170R'), interactive = T)

#' ## P-Value Distribution
dt <- get_total_mae_result(sa, rna_id = '103170R')
dt <- dt[as_gt != '0/0']
hist(dt$pvalue)
qqplot(-log10(ppoints(nrow(dt))), -log10(sort(dt$pvalue)))
abline(0,1)
grid()
hist(dt$alt_freq)
