#'---
#' title: All candidates from Exome, RNA, and Protein data
#' author: Daniel Bader
#' wb:
#'   input: 
#'     - exome_candy: "`sm config['PROC_RESULTS'] + 'variants_wes_candidates.RDS'`"
#'     - rna_aber_exp: "`sm config['PROC_RESULTS'] + 'rna_aberrant_expression.RDS'`"
#'     - rna_splice: "`sm config['PROC_RESULTS'] + 'rna_splicing_fraser_results.RDS'`"
#'     - rna_mae: "`sm config['PROC_RESULTS'] + 'rna_mae_deseq_results.RDS'`"
#'     - prot_aber_exp: "`sm config['PROC_RESULTS'] + 'proteome_aberrant_expression.tsv'`"
#'   output: [
#'     all_candidates: "`sm config['PROC_RESULTS'] + 'all_candidates.tsv'`"
#'   ]
#' output: 
#'   html_document:
#'     toc_float: yes
#'     toc: yes
#'---
#'

#+ echo=F, warning=F, message=F
source("src/r/config.R")
file_all_candidates <- snakemake@output[['all_candidates']]

# very noisy gene anno file
file_disease_gene_anno <- file.path(RAWDIR, "gene_info/meta_disease_genes.tsv")

# tidy results
file_exome_candy <- snakemake@input[['exome_candy']]
file_rna_aber_exp <- snakemake@input[['rna_aber_exp']]
file_rna_splice <- snakemake@input[['rna_splice']]
file_rna_mae <- snakemake@input[['rna_mae']]
file_prot_aber_exp <- snakemake@input[['prot_aber_exp']]


exome_candy <- readRDS(file_exome_candy)
rna_aber_exp <- readRDS(file_rna_aber_exp)
# rna_splice <- readRDS(file_rna_splice)
rna_mae <- readRDS(file_rna_mae)
prot_aber_exp <- fread(file_prot_aber_exp)
GENE_ANNO <- fread(file_disease_gene_anno, na.strings=c('NA',''))


# 
# reduce tables to candidates only
#
gene_anno_diet <- GENE_ANNO[,.(GENE_ID=HGNC_GENE_NAME, DISEASE, MIM_NUMBERS)]
exome_diet_candy <- exome_candy[,.(EXOME_ID=sample, FIBROBLAST_ID, GENE_ID=hgncid, EXOME_IS_SIGNI=TRUE)]
rna_diet_aber_exp <- rna_aber_exp[rna_is_signi==TRUE, .(FIBROBLAST_ID, GENE_ID=hgncid, RNA_IS_ABER_EXP=rna_is_signi)]
rna_diet_mae <- rna_mae[MAE_IS_SIGNI==TRUE, .(FIBROBLAST_ID, GENE_ID=HGNCID, MAE_IS_SIGNI)]
prot_diet_aber_exp <- prot_aber_exp[PROT_IS_ABER_EXP==TRUE, .(FIBROBLAST_ID, GENE_ID=GENE_NAME, PROT_IS_ABER_EXP)]

#
# merge
#
merged_candy_dt <- merge(merge(merge(
    merge(exome_diet_candy, rna_diet_aber_exp, all=T), 
    rna_diet_mae, all=T),
    prot_diet_aber_exp, all=T),
    gene_anno_diet, by='GENE_ID', all.x=T
)

merged_candy_dt <- merge(
    merged_candy_dt, 
    SAMPLE_ANNOTATION[,.(FIBROBLAST_ID, EXOME_ID, DIAGNOSED_GENE=KNOWN_MUTATION)],
    by=c('FIBROBLAST_ID', 'EXOME_ID'),
    all.x=T
)

# SAVE
write_tsv(merged_candy_dt, file = file_all_candidates)



#' ## Results 
#' 
#' - Exome variant prioritization: rare, protein-affecting, potentially bi-allelic
#' - RNA aberrant expression: absolute Z-score > 3, adjusted P-value < 0.05
#' - RNA mono-allelic expression: adjusted P-value < 0.05, alt_allele_freq > 0.8
#' - Protein aberrant expression: absolute Z-score > 3, adjusted P-value < 0.05
#' 
#' Missing: splicing from FraseR, waiting for cutoffs
#' 
#' ## Special columns
#' 
#' - DISEASE = [MITO, MAYBE, OTHER, NA]: annotated disease from HMGU tables, to be extended with more terms from OMIM
#' - MIM_NUMBERS: If annotated the numbers for OMIM entries.
#' - DIAGNOSED_GENE: "NA" for an unsolved patient, or gene name
#' 


## STRONG CANDIDATES
#+
DT::datatable(
    merged_candy_dt[is.na(DIAGNOSED_GENE) & !is.na(DISEASE)], 
    filter='top', 
    rownames = FALSE
    # ,options = list(scrollX = TRUE)
)



## All candidates
#+ echo=F
#head(exome_display_dt)
DT::datatable(
    merged_candy_dt, 
    filter='top', 
    rownames = FALSE
    # ,options = list(scrollX = TRUE)
)


