#'---
#' title: Tidy proteomics from Pichler 
#' author: Daniel Bader
#' wb:
#'   input: /s/project/mitoMultiOmics/raw_data/proteome/20160202_pichler_proteome/20151005_QExacitveHF_100min_proteinGroups.txt
#'   output: "/s/project/patient_report/tidy_results//proteome_pichler_100min.tsv"
#' output: 
#'   html_document:
#'     toc_float: yes
#'     toc: yes
#'---
#'

#+ echo=F
source("src/r/config.R")

#' # FILE 

#+
file_pichler_100min <- file.path(
    RAWDIR, 
    "proteome/20160202_pichler_proteome/",
    "20151005_QExacitveHF_100min_proteinGroups.txt"
)
stopifnot(file.exists(file_pichler_100min))

print(file_pichler_100min)

#' # Read and process data

pdt <- wrapper_proteinGroupsTxt_to_tidy_table(
    file_pichler_100min, 
    intensity_column_pattern = "LFQ.intensity."
)
print(pdt)


#' # Save tidy table

file_tidy_pichler_100min <- file.path(
    TIDYDIR,
    "proteome_pichler_100min.tsv"
)

write_tsv(pdt, file=file_tidy_pichler_100min)

stopifnot(file.exists(file_tidy_pichler_100min))
