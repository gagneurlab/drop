#'---
#' title: RNA candidates by aberrant and mono-allelic expression
#' author: Daniel Bader
#' wb:
#'   input: 
#'     - rna_aber_exp: "`sm config['PROC_RESULTS'] + 'rna_aberrant_expression.RDS'`"
#'     - rna_mae: "`sm config['PROC_RESULTS'] + 'rna_mae_deseq_results.RDS'`"
#'   output: 
#' output: 
#'   html_document
#'---
#'

#+ echo=F
source("src/r/config.R")
file_rna_aber_exp <- snakemake@input[['rna_aber_exp']]
file_rna_mae <- snakemake@input[['rna_mae']]
file_disease_gene_anno <- file.path(RAWDIR, "gene_info/meta_disease_genes.tsv")

GENE_ANNO <- fread(file_disease_gene_anno, na.strings=c('NA',''))
rna_aber_exp <- readRDS(file_rna_aber_exp)
rna_mae <- readRDS(file_rna_mae)


#' 
#' ## RNA aberrant expression
#' 

#+ echo=F
# subset columns
columns_to_show_aberexp <- c(
    "FIBROBLAST_ID", 
    "hgncid",
    "log2FoldChange",
    "rna_zscore",
    "rna_hochberg_padj",
    "rna_norm_lvl",
    "rna_is_signi"
)


# merge exome and disease gene info
display_dt <- merge(
    GENE_ANNO[,.(HGNC_GENE_NAME, MIM_NUMBERS, DISEASE)], 
    rna_aber_exp[, ..columns_to_show_aberexp], 
    by.x='HGNC_GENE_NAME', 
    by.y='hgncid', 
    all.y=T
)

# round columns with numbers
columns_signif <- c('log2FoldChange', "rna_zscore", "rna_hochberg_padj", 'rna_norm_lvl')
for(j in columns_signif){
    display_dt[, c(j):= list(signif(get(j), digits = 3))]
}


#'  
#'  RNA aberrant expression: absolute Z-score > 3, adjusted P-value < 0.05
#' 
#+ echo=F
#head(display_dt)
DT::datatable(
    display_dt, 
    filter='top', 
    rownames = FALSE
    # ,options = list(scrollX = TRUE)
)


#' 
#' ## RNA mono-allelic expression
#' 

#+ echo=F
# subset columns
columns_to_show_mae <- setdiff(names(rna_mae), 'pvalue')


# merge exome and disease gene info
display_dt <- merge(
    GENE_ANNO[,.(HGNC_GENE_NAME, MIM_NUMBERS, DISEASE)], 
    rna_mae[, ..columns_to_show_mae],
    by.x='HGNC_GENE_NAME', 
    by.y='HGNCID', 
    all.y=T
)

# round columns with numbers
columns_signif <- c('padj', "exacmaf", "alt_allele_freq")
for(j in columns_signif){
    display_dt[, c(j):= list(signif(get(j), digits = 3))]
}



#' 
#' RNA mono-allelic expression: adjusted P-value < 0.05, alt_allele_freq > 0.8
#' 
#+ echo=F
#head(display_dt)
DT::datatable(
    display_dt, 
    filter='top', 
    rownames = FALSE
    # ,options = list(scrollX = TRUE)
)



#+ END, echo=F


