#'---
#' title: Mito Disease Score
#' author: vyepez
#' wb:
#'  input:
#'   - ods_ss: '`sm config["PROC_RESULTS"] + "/v29_overlap/outrider/fib_ss/ods.Rds"`'
#'   - ods_ns: '`sm config["PROC_RESULTS"] + "/v29_overlap/outrider/fib_ns/ods.Rds"`'
#'   - exomes_table: "../exomes1000/processed_data/exomes_clean.tsv"
#'  output:
#' output:
#'   html_document:
#'    code_folding: show
#'    code_download: TRUE
#'---

#+ echo=F
saveRDS(snakemake, "tmp/mdc_score.snakemake")
# snakemake <- readRDS("tmp/mdc_score.snakemake")

suppressPackageStartupMessages({
    library(OUTRIDER)
    library(ggplot2)
    library(data.table)
    library(dplyr)
    library(DESeq2)
})

source("Scripts/_functions/gene_annotation/add_gene_info_cols.R")
source("../exomes1000/Scripts/_functions/useful_functions.R")  # needed for duplicate_row()

#' ## Read, normalize and center data
ods_ss <- readRDS(snakemake@input$ods_ss)
ods_ns <- readRDS(snakemake@input$ods_ns)

norm_sizefactor_log <- function(ods){
    counts <- counts(ods, normalized = F)
    counts <- t(t(counts) / estimateSizeFactorsForMatrix(counts))
    counts <- log2(counts + 1) - rowMeans2(log2(counts + 1))
    return(counts)
}

counts_ss <- norm_sizefactor_log(ods_ss)
counts_ns <- norm_sizefactor_log(ods_ns)

# Read exo_dt
exo_dt <- fread(snakemake@input$exomes_table)
exo_dt <- exo_dt[WES_in_db == 'Y', .(RNA_ID, MDC_score)] %>% na.omit()
exo_dt <- exo_dt[RNA_ID != '.']
# Some RNAs have either ID_1 & ID_2, or ID_1, ID_2
exo_dt <- duplicate_row(exo_dt, sep = '&', col = 'RNA_ID')
exo_dt <- duplicate_row(exo_dt, sep = ',', col = 'RNA_ID')
cols <- names(exo_dt)
exo_dt[, (cols) := lapply(.SD, function(x) gsub(" ", "", x))]  # Remove spaces
exo_dt[, MDC_score := as.numeric(MDC_score)]

# exo_dt <- exo_dt[MDC_score < 8]

#' Samples in 1000exomes, but not in RNA
setdiff(exo_dt$RNA_ID, c(colnames(ods_ss), colnames(ods_ns))) %>% sort
#' Samples in RNA, but not in 1000exomes
setdiff(c(colnames(ods_ss), colnames(ods_ns)), exo_dt$RNA_ID) %>% sort

#' ## Compute the correlation
#' Number of samples with both MDC score and RNA
get_cor_exp_mds_dt <- function(COUNTS, EXO_DT){
    
    inter_samples <- sort(intersect(colnames(COUNTS), EXO_DT$RNA_ID))
    print(inter_samples)
    
    COUNTS <- COUNTS[, inter_samples]
    EXO_DT <- EXO_DT[RNA_ID %in% inter_samples]
    setorder(EXO_DT, RNA_ID)
    
    cort <- apply(COUNTS, 1, 
              function(x){
                  ct <- cor.test(x, EXO_DT[, MDC_score], method = 'spearman')
                  return(c(pv = ct$p.value, ct$estimate))
                  }
              )
    cor_dt <- t(cort) %>% as.data.table(keep.rownames = T)
    setnames(cor_dt, 'rn', 'gene')
    cor_dt[, padj := p.adjust(pv, method = 'BH')]
    cor_dt <- add_hans_class(cor_dt, gene_name_col = 'gene', return_all_info = F)
    setorder(cor_dt, padj)
    return(cor_dt)
}


ct_ss <- get_cor_exp_mds_dt(counts_ss, exo_dt)
ct_ns <- get_cor_exp_mds_dt(counts_ns, exo_dt)
ct_ns[, MITO_DISEASE_GENE := NULL]
ct_ss[padj < .5]
ct_ns[padj < .5]

# Combine tables and add fisher pv
ct_all <- inner_join(ct_ss, ct_ns, by = 'gene') %>% as.data.table
library(metap)
ct_all$fisher_pv = sapply(1:nrow(ct_all), function(i) sumlog(unlist(ct_all[i, .(pv.x, pv.y)]))$p)  # sumlog is not vectorized
ct_all[sign(rho.x) != sign(rho.y), fisher_pv := 1]   # Fisher pvalue valid only if the correlations had the same sign
setorder(ct_all, fisher_pv)
ct_all = add_mitocarta_col(ct_all, gene_name_col = 'gene')

DT::datatable(ct_all[fisher_pv < .01], style = 'bootstrap')

#' ## Visualize top hits
plot_cor_rna_mds <- function(gene, COUNTS_SS, COUNTS_NS, EXO_DT){
    inter_samples <- sort(intersect(colnames(COUNTS_SS), EXO_DT$RNA_ID))
    COUNTS_SS <- COUNTS_SS[, inter_samples]
    EXO <- EXO_DT[RNA_ID %in% inter_samples]
    setorder(EXO, RNA_ID)
    plot_ss = data.table(counts = COUNTS_SS[gene, ], mdc = EXO$MDC_score, STRANDED = 'Specific')
    
    inter_samples <- sort(intersect(colnames(COUNTS_NS), EXO_DT$RNA_ID))
    COUNTS_NS <- COUNTS_NS[, inter_samples]
    EXO <- EXO_DT[RNA_ID %in% inter_samples]
    setorder(EXO, RNA_ID)
    plot_ns = data.table(counts = COUNTS_NS[gene, ], mdc = EXO$MDC_score, STRANDED = 'Non Specific')
    
    plot_dt <- rbind(plot_ss, plot_ns)
    
    g <- ggplot(plot_dt, aes(as.factor(mdc), counts)) + geom_boxplot() + geom_beeswarm() + facet_wrap(~ STRANDED) + 
        geom_smooth(aes(eval(mdc+1), counts), method = 'lm') + 
        labs(title = gene, x = 'Mito Disease Criteria Score', y = 'L2FC') + theme_bw()
    
    return(g)
}

library(ggbeeswarm)
plot_cor_rna_mds('MIEF1', counts_ss, counts_ns, exo_dt)
plot_cor_rna_mds('RGPD2', counts_ss, counts_ns, exo_dt)
plot_cor_rna_mds('TMEM164', counts_ss, counts_ns, exo_dt)
plot_cor_rna_mds('SLC25A44', counts_ss, counts_ns, exo_dt)
plot_cor_rna_mds('BID', counts_ss, counts_ns, exo_dt)



