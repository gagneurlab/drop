#'---
#' title: Disease associated genes
#' author: Daniel Bader, Vicente Yépez
#' wb:
#'   input: "/s/project/mitoMultiOmics/raw_data/gene_info//meta_disease_genes.tsv"
#'   output: 
#' output: 
#'   html_document
#'---
#'

#+ echo=F
library(data.table)
library(DT)

dir_gene_info <- '/s/project/mitoMultiOmics/raw_data/gene_info/'

#' Hans list of genes
cleaned_file <- file.path(dir_gene_info, 'mitochondrial_disorder_genes_prokisch_mayr_cleaned.tsv')
prokisch_mayr_dt <- fread(cleaned_file)

# New order and subset
coln_order <- c(
    'HGNC_GENE_NAME', 
    'DISEASE', 
    'CATEGORY',
    'MIM_NUMBERS', 
    'YEAR',
    "CARDIOMYOPATHY",
    "FUNCTION",
    "CELLULAR_LOCALISATION",
    "ASSOCIATED_DISEASE_PHENOTYPES",
    "OTHER_DESIGNATIONS",
    'LOCUS', 
    'ENTREZ_GENE_ID', 
    'OTHER_ALIASES'
)


prokisch_mayr_dt <- prokisch_mayr_dt[, coln_order, with=F]
DT::datatable(prokisch_mayr_dt, filter='top', rownames = FALSE)

barplot(table(prokisch_mayr_dt$DISEASE), las = 1)

#' Mito disease genes table
file_disease_gene <- file.path(dir_gene_info, 'disease_genes_from_sample_anno.tsv')
genetable <- fread(file_disease_gene)

DT::datatable(genetable, filter='top', rownames = FALSE)
barplot(table(genetable$DISEASE), las = 2)



