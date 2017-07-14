#'---
#' title: R script
#' author: Daniel Bader
#' wb:
#'   input: [
#'     "/s/project/patient_report/tidy_results/proteome_pichler_100min.tsv"
#'   ]
#'   output: "/s/project/patient_report/tidy_results/roteome_aberrant_expression.tsv"
#' output: 
#'   html_document:
#'     toc_float: yes
#'     toc: yes
#'---
#'



file_tidy_pichler_100min <- file.path(
    TIDYDIR,
    "proteome_pichler_100min.tsv"
)
file_aber_prot_exp <- file.path(
    TIDYDIR,
    "proteome_aberrant_expression.tsv"
)

pdt <- fread(file_tidy_pichler_100min)


#'
#' # Aberrant expression
#'

prot_intensity <- convert_tidy_table_to_numeric_matrix(
    pdt, 'GENE_NAME', 'PROTEOME_ID', 'LFQ_INTENSITY'
)
# head(prot_intensity)

#' 
#' ## Check clustering
#' 
# prot_log2fc_intensity = normalize_expression_matrix(
#    prot_intensity, sizefactor = T, rowcenter = T, log2scale = T, nonzero = 0
# )
# d <- dist(prot_log2fc_intensity)
# hc <- hclust(d)
# heatmap_notrace()


#' 
#' ## Compute DE with limma
#' 
#' * transform protein intensity to log2 space
#' * normalize with DESeq size factors
#' * generate design matrix
#' * check if some samples are not respected in the design
#' * limma DE fit
#' * compute Z-score
#' 
prot_aberexp <- wrapper_aberrant_protein_expr(prot_intensity)
head(prot_aberexp)

#'
#' ## SAVE
#'
file_aber_prot_exp
write_tsv(prot_aberexp, file = file_aber_prot_exp)

