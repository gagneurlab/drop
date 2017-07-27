#'---
#' title: Disease associated genes
#' author: Daniel Bader
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
file_meta_disease_gene <- file.path(dir_gene_info, 'meta_disease_genes.tsv')
genetable <- fread(file_meta_disease_gene)

coln_order <- c(
    'HGNC_GENE_NAME', 
    'DISEASE', 
    'MIM_NUMBERS', 
    'OMIM_LINK', 
    'YEAR',
    "SUBGROUP",
    "NEURO",
    "LEIGH",
    "FUNCTION",
    "CELLULAR_LOCALISATION",
    "AFFECTED_PATHWAY",
    "ANNOTATION",
    "ASSOCIATED_DISEASE_PHENOTYPE_S",
    "OTHER_DESIGNATIONS",
    'LOCUS', 
    'ENTREZ_GENE_ID', 
    'OTHER_ALIASES' 
    )
genetable <- genetable[, coln_order, with=F]


DT::datatable(genetable, filter='top', rownames = FALSE)


