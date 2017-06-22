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

#protdir <- file.path(RAWDIR, "proteome", "20170614_kopajtich_kuester_proteome")
protdir <- file.path("~/Documents/tmp_kuester_proteome")

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

#' Create tidy proteome data table "pdt"
pdt <- read_data_table_from_proteinGroupsTxt(
    tmp_file_kuester, 
    intensity_column_pattern = "LFQ.intensity.", 
    column_ids = c("Protein.IDs", "Gene.names"), 
    column_sample_id = "PROTEOME_ID", 
    column_intensity = "LFQ_INTENSITY"
)
print(head(pdt))


#' Filter protein expression by gene properties
#' 
pdt <- filter_proteome_by_gene_properties(
    pdt, verbose=T
)
print(pdt)


#' Compute NA frequencies
#' 
compute_na_frequency(pdt, column_id = 'GENE_NAME')
compute_na_frequency(pdt, column_id = 'PROTEOME_ID')
print(pdt)


#' 
#' ## iBAQ
#'

#' Extract intensities from proteinGroups.txt
protein_table_ibaq <- proteinGroups_read_table( 
    tmp_file_kuester, 
    intensity_column_pattern = "iBAQ"
)
protein_table_ibaq



