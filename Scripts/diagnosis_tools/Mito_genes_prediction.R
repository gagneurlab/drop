#'---
#' title: Mito genes prediction
#' author: vyepez
#' wb:
#'  input:
#'   - pred_scores: '../mito-ncRNA/results/gtex_leostable_OCR.csv'
#'  output:
#' output:
#'   html_document:
#'    code_folding: show
#'    code_download: TRUE
#'---

#+ echo=F
saveRDS(snakemake, "tmp/mito_pred_scores.snakemake")
# snakemake <- readRDS("tmp/mito_pred_scores.snakemake")

suppressPackageStartupMessages({
    library(ggplot2)
    library(ggthemes)
    library(data.table)
    library(magrittr)
    library(cowplot)
    library(ggpval)
})

#' ## Read the prediction table
pred_dt <- fread(snakemake@input$pred_scores)
#' Total number of genes
nrow(pred_dt)

#' Total number of protein coding genes that are not MITOCARTA
pred_dt <- pred_dt[y == FALSE]
nrow(pred_dt)

setnames(pred_dt, c("y", "prob.y"), c("MITOCARTA", "prob_mito"))

setorder(pred_dt, -prob_mito)

#' ## Plot the distribution
ggplot(pred_dt, aes(prob_mito)) + geom_histogram() + theme_cowplot() + ggtitle("Mito score distribution")

#' Add genes diagnosed from 2000exomes
exo_db <- fread("../exomes1000/processed_data/exomes_clean.tsv")
exo_db <- exo_db[WES_in_db == 'Y']

pred_dt[, in2000exomes := gene_name %in% exo_db$gene_name]
table(pred_dt$in2000exomes)
table(pred_dt[prob_mito > 0.2, in2000exomes])
pred_dt[prob_mito > 0.2 & in2000exomes == TRUE]

g <- ggplot(pred_dt, aes(in2000exomes, prob_mito)) + geom_boxplot() + theme_cowplot()
ggpval::add_pval(g, pairs = list(c(1,2)), heights = .8, textsize = 4)

#' Check if genes are expressed in our fibroblasts
ods_ss <- readRDS("/s/project/genetic_diagnosis/processed_results/v29_overlap/outrider/ss/ods.Rds")
pred_dt[, expressed := gene_name %in% toupper(row.names(ods_ss))]
table(pred_dt$expressed)
g <- ggplot(pred_dt, aes(expressed, prob_mito)) + geom_boxplot() + theme_cowplot()
ggpval::add_pval(g, pairs = list(c(1,2)), heights = .8, textsize = 4)

#' ## Download the table
#+ echo=F
DT::datatable(pred_dt)
write.table(pred_dt, "/s/public_webshare/project/genetic_diagnosis/results/mito_pred_scores.txt", sep = "\t", quote = F, row.names = F)
#' [Download prediction scores table](https://i12g-gagneurweb.informatik.tu-muenchen.de/project/genetic_diagnosis/results/mito_pred_scores.txt)


