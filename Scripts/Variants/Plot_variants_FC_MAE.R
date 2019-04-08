#'---
#' title: Variants FC and MAE plots
#' author: vyepez
#' wb:
#'  input:
#'   - variants_outrider: '`sm config["PROC_RESULTS"] + "/process_vcf/vars_outrider/outrider_all_vars.Rds"`'
#'   - mae_outrider: '`sm config["PROC_RESULTS"] + "/processed_results/mae/mae_results_outrider.Rds"`'
#'  output:
#'  threads: 10
#' output:
#'   html_document:
#'    code_folding: show
#'    code_download: TRUE
#'---

#+echo=F
saveRDS(snakemake, "tmp/outrider_vars_plot.snakemake")
# snakemake <- readRDS("tmp/outrider_vars_plot.snakemake")

suppressPackageStartupMessages({
  library(data.table)
  library(magrittr)
  library(ggplot2)
  library(cowplot)
  library(ggthemes)
  library(ggbeeswarm)
})


#' ## Gene Expression FC vs variant class
pt <- readRDS(snakemake@input$variants_outrider) 

# Violin plot
#+ fig.width=8
ggplot(pt[!is.na(AF_CAT)], aes(var_type, FC)) + geom_violin(aes(fill = AF_CAT)) + geom_hline(yintercept = 1) + 
  scale_y_log10(limits = c(.1, 10)) + scale_fill_ptol() + facet_wrap(~ genotype, nrow = 2)

# Boxplot
#+ fig.width=8
ggplot(pt[!is.na(AF_CAT)], aes(var_type, FC)) + geom_boxplot(aes(fill = AF_CAT)) + geom_hline(yintercept = 1) + 
  scale_y_log10(limits = c(.5, 2)) + scale_fill_ptol() + facet_wrap(~ genotype, nrow = 2)

# Boxplot, protein coding genes only
#+ fig.width=8
ggplot(pt[!is.na(AF_CAT) &  gene_type == 'protein_coding'], aes(var_type, FC)) + geom_boxplot(aes(fill = AF_CAT)) + geom_hline(yintercept = 1) + 
  scale_y_log10(limits = c(.5, 2)) + scale_fill_ptol() + facet_wrap(~ genotype, nrow = 2)


# outliers only
#+ fig.width=8
ggplot(pt[!is.na(AF_CAT) & aberrant == T], aes(var_type, FC, color = AF_CAT)) + # geom_violin(aes(fill = AF_CAT)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_quasirandom(dodge.width = .8) + 
  geom_hline(yintercept = 1) + 
  scale_y_log10() + 
  scale_color_ptol() + facet_wrap(~ genotype, nrow = 2) + ggtitle("Outliers only")



#' ## MAE
mae_res_all <- readRDS(snakemake@input$mae_outrider)

# Violin
#+ fig.width=8
ggplot(mae_res_all[!is.na(AF_CAT)], aes(var_type, allele_dif)) + geom_violin(aes(fill = AF_CAT)) + 
  scale_fill_ptol() + geom_hline(yintercept = median(mae_res_all$allele_dif)) + labs(y = 'ref freq - alt freq')

# Boxplot
#+ fig.width=8
ggplot(mae_res_all[!is.na(AF_CAT)], aes(var_type, allele_dif)) + geom_boxplot(aes(fill = AF_CAT)) + 
  scale_fill_ptol() + geom_hline(yintercept = median(mae_res_all$allele_dif), color = 'gray') + labs(y = 'ref freq - alt freq')


# Outliers only
#+ fig.width=8
ggplot(mae_res_all[!is.na(AF_CAT) & aberrant == T], aes(var_type, allele_dif, color = AF_CAT)) + 
  geom_boxplot(aes(fill = AF_CAT), outlier.shape = NA) + 
  geom_quasirandom(dodge.width = .8) + 
  geom_hline(yintercept = median(mae_res_all$allele_dif), color = 'gray') + 
  scale_color_grey() + scale_fill_brewer(palette = 'Blues', direction = -1) + 
  labs(title = "Outliers only", y = 'ref freq - alt freq')

# Check both samples with rare stop mutations and alternative higher
DT::datatable(mae_res_all[!is.na(AF_CAT) & var_type == 'stop' & allele_dif < 0 & aberrant == T & MAX_AF < .001])
