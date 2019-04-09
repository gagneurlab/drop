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

#' ## (Rare) variants associated with outliers
pt <- readRDS(snakemake@input$variants_outrider) 

# Add outlier class
pt[, outlier_class := 'non-outlier']
pt[aberrant == T & FC > 1, outlier_class := 'overexpression']
pt[aberrant == T & FC < 1, outlier_class := 'underexpression']

# Plot all variants
pd <- pt[, .N, by = .(outlier_class, var_type)]
pd[, prop := N/sum(N), by = outlier_class]
totals <- pd[, sum(N), by = outlier_class]
pd[, var_type := factor(var_type, levels = c("no rare variant", "non_coding", "synonymous", "coding", "frameshift", "splice", "stop"))]

#+ fig.width=10, fig.height=3
ggplot(pd, aes(outlier_class, prop)) + geom_bar(stat= 'identity', aes(fill = var_type)) + scale_fill_ptol() +
  coord_flip() +
  scale_y_continuous(labels=scales::percent) +
  geom_text(aes(outlier_class, 1.05, label = V1, fill = NULL), data = totals) +
  labs(y = "", title = "All variants")

# Plot rare variants only
pt[, rare_var_type := var_type]
pt[MAX_AF > .01, rare_var_type := "no rare variant"]
pr <- pt[, .N, by = .(outlier_class, rare_var_type)]
pr[, prop := N/sum(N), by = outlier_class]
setorder(pr, outlier_class)
totals <- pr[, sum(N), by = outlier_class]
pr[, rare_var_type := factor(rare_var_type, levels = c("no rare variant", "non_coding", "synonymous", "coding", "frameshift", "splice", "stop"))]

#+ fig.width=10, fig.height=3
ggplot(pr, aes(outlier_class, prop)) + geom_bar(stat= 'identity', aes(fill = rare_var_type)) + scale_fill_ptol() +
  coord_flip() +
  scale_y_continuous(labels=scales::percent) +
  geom_text(aes(outlier_class, 1.05, label = V1, fill = NULL), data = totals) +
  labs(y = "", title = "Rare variants only")


#' ## Gene Expression FC vs variant class

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
