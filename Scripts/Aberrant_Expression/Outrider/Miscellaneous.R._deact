#'---
#' title: Miscellaneous plots
#' author: vyepez
#' wb:
#'  input:
#'   - ods_ss: '`sm expand(config["PROC_RESULTS"] + "/{annotation}/outrider/fib_ss/ods.Rds", annotation=config["ANNOTATIONS"])`'
#'   - ods_ns: '`sm expand(config["PROC_RESULTS"] + "/{annotation}/outrider/fib_ns/ods.Rds", annotation=config["ANNOTATIONS"])`'
#'  output:
#' output:
#'   html_document:
#'    code_folding: show
#'    code_download: TRUE
#'---

#+ echo=F
saveRDS(snakemake, "tmp/outrider_misc.snakemake")
# snakemake <- readRDS("tmp/outrider_misc.snakemake")

suppressPackageStartupMessages({
    library(OUTRIDER)
    library(ggplot2)
    library(data.table)
    library(dplyr)
    library(ggbeeswarm)
})


#'
ods_ss <- readRDS(snakemake@input$ods_ss)
dim(ods_ss)

#' Mirjana is interested in seeing the expression of the following mito genes on sample 99590R
mito_genes <- c("MT-CO1","MT-CO2","MT-CYB","MT-RNR2","MT-ND1","MT-ND2","MT-ND3","MT-ND5","MT-ND6")
other_genes <- c("MTERF1", "NRF1", "NFE2L2") # MTERF1 should be down; PRKAA2 didn't pass the filter
sample = '99590R'

m_raw <- melt(counts(ods_ss[c(mito_genes, other_genes),], normalized = F)) %>% as.data.table
colnames(m_raw) <- c("gene", "RNA_ID", "raw_counts")

m_norm <- melt(counts(ods_ss[c(mito_genes, other_genes),], normalized = T)) %>% as.data.table
colnames(m_norm) <- c("gene", "RNA_ID", "norm_counts")

mt_genes_dt <- cbind(m_raw, m_norm[, .(norm_counts)])
mt_genes_dt[, size_factors := rep(sizeFactors(ods_ss), each = length(c(mito_genes, other_genes)))]
mt_genes_dt[, sF_counts := raw_counts / size_factors]
mt_genes_dt[, int_sample := F]
mt_genes_dt[RNA_ID == sample, int_sample := T]

mt_sub <- mt_genes_dt[gene %in% c("MT-ND1", "MT-ND5", "MT-ND6", "MTERF1", "NRF1")]
mt_sub[, mt_genes_dt := factor(mt_genes_dt, levels = c("MT-ND1", "MT-ND5", "MT-ND6", "MTERF1", "NRF1"))]

#' sF_counts: counts normalized by sequencing depth using size factors

# beanplots
expression_beanplot <- function(DT, counts_type = c('raw_counts', 'norm_counts', 'sF_counts'), main = ''){
    ggplot(DT, aes(gene, get(counts_type))) + geom_beeswarm(aes(col = int_sample)) + 
    scale_y_log10() + theme_bw(base_size = 14) + theme(legend.position="none") +
    scale_color_manual(values = c("gray60", "firebrick")) + labs(y = counts_type, title = main)
}
#+ fig.width=11
expression_beanplot(mt_genes_dt[gene %in% mito_genes], counts_type = 'raw_counts')
expression_beanplot(mt_genes_dt[gene %in% mito_genes], counts_type = 'norm_counts')
#+ fig.width=5
expression_beanplot(mt_genes_dt[gene %in% other_genes], counts_type = 'raw_counts')
expression_beanplot(mt_genes_dt[gene %in% other_genes], counts_type = 'norm_counts')

# boxplots
expression_boxplot <- function(DT, counts_type = c('raw_counts', 'norm_counts', 'sF_counts'), main = ''){
    ggplot(DT, aes(gene, get(counts_type))) + geom_boxplot() + 
    geom_point(data = DT[RNA_ID == sample], aes(gene, get(counts_type)), color = 'firebrick') + 
    scale_y_log10() + theme_bw(base_size = 14) + labs(y = counts_type, title = main)
}

#+ fig.width=11
expression_boxplot(mt_genes_dt[gene %in% mito_genes], counts_type = 'raw_counts')
expression_boxplot(mt_genes_dt[gene %in% mito_genes], counts_type = 'sF_counts')
expression_boxplot(mt_genes_dt[gene %in% mito_genes], counts_type = 'norm_counts')
#+ fig.width=5
expression_boxplot(mt_genes_dt[gene %in% other_genes], counts_type = 'raw_counts')
expression_boxplot(mt_genes_dt[gene %in% other_genes], counts_type = 'sF_counts')
expression_boxplot(mt_genes_dt[gene %in% other_genes], counts_type = 'norm_counts')

#+ fig.width=6
expression_boxplot(mt_sub, counts_type = 'sF_counts')

#' Volcano plot of sample
plotVolcano(ods_ss, sample)

#'
x = counts(ods_ss)
library(LSD)
heatscatter(x["SMIM26", ], x["SURF1", ], log = 'xy', cor = T, xlab = 'SMIM26', ylab = 'SURF1'); grid(); abline(0,1)
