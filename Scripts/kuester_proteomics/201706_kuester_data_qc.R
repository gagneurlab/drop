#'---
#' title: Kuester data QC
#' author: Daniel Bader
#' wb:
#'   input: 
#'   output: 
#' output: 
#'   html_document:
#'     toc_float: yes
#'     toc: yes
#'---
#'

#+ echo=F
source("src/r/config.R")

#' # Data

protdir <- file.path(RAWDIR, "proteome", "20170614_kopajtich_kuester_proteome")
files_kuester <- list.files(protdir, pattern = "^m.*txt$", full.names = T)

#' Kuester proteome files:
print(files_kuester)



#' 
#' # Single Shot
#' 
tmp_file_kuester <- files_kuester[1]

#' 
#' ## LFQ
#'


#' Extract intensities from proteinGroups.txt
protein_table <- proteinGroups_read_table( 
    tmp_file_kuester, 
    intensity_column_pattern = "LFQ" 
)
print(head(protein_table))

#' Create expression matrix
#+ message='show'
mat <- proteinGroups_quality_control(
    protein_table, 
    intensity_column_pattern = "LFQ", 
    max_na_frequency = RELIABLE_PROT_FRACTION_NA, 
    low_expr_quantile = LOW_EXPR_QUANTILE,
    return_matrix = TRUE
)
head(mat)


#' 
#' ## iBAQ
#'




