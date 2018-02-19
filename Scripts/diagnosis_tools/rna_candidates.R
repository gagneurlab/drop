#'---
#' title: RNA candidates by aberrant and mono-allelic expression
#' author: Daniel Bader, Vicente Yépez
#' wb:
#'   input: 
#'     - rna_aber_exp: "`sm config['PROC_RESULTS'] + 'rna_aberrant_expression.RDS'`"
#'     - rna_aber_b1_b2: "`sm config['PROC_RESULTS'] + 'rna_aberrant_expression_b1_b2.RDS'`"
#'     - rna_aber_b3: "`sm config['PROC_RESULTS'] + 'rna_aberrant_expression_b3.RDS'`"
#'     - rna_mae: "`sm config['PROC_RESULTS'] + 'rna_mae_deseq_results.RDS'`"
#'   output: 
#' output: 
#'   html_document
#'---
#'


#+ echo=F
source("src/r/config.R")
file_rna_aber_exp <- snakemake@input[['rna_aber_exp']]
file_rna_aber_b1_b2 <- snakemake@input[['rna_aber_b1_b2']]
file_rna_aber_b3 <- snakemake@input[['rna_aber_b3']]
file_rna_mae <- snakemake@input[['rna_mae']]
file_disease_gene_anno <- file.path(RAWDIR, "gene_info/meta_disease_genes.tsv")

# file_rna_aber_exp = file.path(PROC_RESULTS, "rna_aberrant_expression.RDS")

GENE_ANNO <- fread(file_disease_gene_anno, na.strings=c('NA',''))
rna_aber_exp <- readRDS(file_rna_aber_exp)
rna_mae <- readRDS(file_rna_mae)

rna_aber_b1_b2 = readRDS(file_rna_aber_b1_b2)
rna_aber_b3 = readRDS(file_rna_aber_b1_b2)
b_all <- rbind(rna_aber_b1_b2, rna_aber_b3) %>% as.data.table()


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

display_dt <- display_dt[rna_is_signi == T]

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


# Batches 1, 2 and 3
# Add OMIM data
omim_dt <- readRDS("/s/project/genetic_diagnosis/resource/omim_dt.Rds")
b_all[, c("baseMean", "lfcSE", "stat", "rawCounts", "log2FoldChange") := NULL]
b_all[, aberrant := abs(zScore) > 3 & padj < .05]
b_all = b_all[aberrant == TRUE]
b_all[, N := .N, by = sampleId]
b_all[, type := ifelse(zScore < 0 , "low", "high")]
disp_new_batches <- left_join(b_all, unique(omim_dt[SYMBOL != "", .(SYMBOL, GMIM, PINH, PMIM)]), 
                  by = c("geneId" = "SYMBOL")) %>% as.data.table()
duplicated(disp_new_batches[, .(sampleId, geneId)]) %>% sum  # Check if there are duplicated due to MIM numbers
disp_new_batches <- disp_new_batches[! (geneId == "CD99" & GMIM == 450000)]
disp_new_batches[, GMIM := as.character(GMIM)]
disp_new_batches[geneId == "CD99", GMIM := "313470, 450000"]

# Add sample annotation data
b1 <- fread("../../LRZ Sync+Share/LMU/sample_annotation/batch1.txt")
b2 <- fread("../../LRZ Sync+Share/LMU/sample_annotation/batch2.txt")
b3 <- fread("../../LRZ Sync+Share/LMU/sample_annotation/batch3.txt")

b123 <- rbind(b1, b2, b3, fill = T)

disp_new_batches <- left_join(disp_new_batches, b123[, .(RNA_ID, CANDIDATE_GENE, 
                               BIOCHEMICAL_DEFECT, CLINICAL_SYMPTOMS, RNA_PERSON)],
                  by = c("sampleId" = "RNA_ID")) %>% as.data.table

setnames(disp_new_batches, old = "sampleId", "RNA_ID")

# round columns with numbers
columns_signif <- c("zScore", "pvalue", "padj")
for(j in columns_signif){
    disp_new_batches[, c(j):= list(signif(get(j), digits = 3))]
}
disp_new_batches[, normalizedCounts := round(normalizedCounts, digits = 2)]

#'  
#'  RNA aberrant expression: absolute Z-score > 3, adjusted P-value < 0.05, New Batches
#' 
#+ echo=F
#head(display_dt)
DT::datatable(
    disp_new_batches, 
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


