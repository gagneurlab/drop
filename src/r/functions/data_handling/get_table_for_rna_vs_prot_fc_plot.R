#'---
#' title: R script
#' author: Daniel Bader
#' output: 
#'   html_document:
#'     toc_float: yes
#'     toc: yes
#'---
#'



get_table_for_rna_vs_prot_fc_plot <- function(
    rna_dt, 
    proteome_dt,
    rna_gene_id='hgncid',
    rna_sample_id='FIBROBLAST_ID',
    rna_fold_change='log2FoldChange',
    prot_gene_id='hgncid',
    prot_sample_id='FIBROBLAST_ID',
    prot_fold_change='prot_log2fc',
    prot_gene_na_freq='prot_freq_na'
){
    merge(
        rna_dt[,c(rna_sample_id, rna_gene_id, rna_fold_change), with=F],
        proteome_dt[, c(prot_sample_id, prot_gene_id, prot_fold_change, prot_gene_na_freq), with=F],
        by.x=c(rna_sample_id, rna_gene_id),
        by.y=c(prot_sample_id, prot_gene_id)
    )
}






