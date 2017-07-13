#'---
#' title: Tidy proteomics from Pichler 
#' author: Daniel Bader
#' wb:
#'   input: 
#'   output: [
#'     "/s/project/patient_report/tidy_results/proteome_pichler_100min_unfiltered_all_samples.tsv",
#'     "/s/project/patient_report/tidy_results/proteome_pichler_100min.tsv"
#'   ]
#' output: 
#'   html_document:
#'     toc_float: yes
#'     toc: yes
#'---
#'

#+ echo=F
source("src/r/config.R")

# protdir <- file.path(
#     RAWDIR, 
#     "proteome/20160202_pichler_proteome/"
# )
protdir <- file.path(
    BADERDIR
)

file_reused_samples <- file.path(
    RAWDIR, "proteome", "protein_ids_measured_2nd_time_in_pichler100min.txt"
)

file_samples_bad_exome2rna <- file.path(
    RAWDIR, "proteome", "protein_ids_with_bad_exome_to_rna_identity.txt"
)
file_tidy_pichler_100min_raw <- file.path(
    TIDYDIR,
    "proteome_pichler_100min_unfiltered_all_samples.tsv"
)
file_tidy_pichler_100min <- file.path(
    TIDYDIR,
    "proteome_pichler_100min.tsv"
)


protein_ids_reused <- readLines(file_reused_samples)
protein_ids_bad_ex2rna <- readLines(file_samples_bad_exome2rna)


#' 
#' # Pichler QExactiveHF 100min
#' 

#+
file_pichler_100min <- file.path(
    protdir,
    "pichler100min_combined_corrected_proteinGroups.txt"
)
stopifnot(file.exists(file_pichler_100min))
print(file_pichler_100min)


#' 
#' ## Read and process data
#' 
pdt <- wrapper_proteinGroupsTxt_to_tidy_table(
    file_pichler_100min, 
    intensity_column_pattern = "LFQ.intensity."
)
print(pdt)


#' 
#' ## Save tidy table
#' 
write_tsv(pdt, file=file_tidy_pichler_100min_raw)
stopifnot(file.exists(file_tidy_pichler_100min_raw))

#' Samples measured
unique(pdt$PROTEOME_ID)


#' 
#' ## Filter low quality or non-standard samples
#'

#' * Samples measured twice from Garwin Pichler 
#'   with bad quality at second time (pichler100min)
protein_ids_reused


#' * Samples with bad exome to RNA identity
protein_ids_bad_ex2rna

#' * Samples with transduced genes or grown 2 weeks
#' 
protein_ids_non_standard <- c( 
    grep("2w", unique(pdt$PROTEOME_ID), value=T),
    grep("\\.T\\.", unique(pdt$PROTEOME_ID), value=T)
)
protein_ids_non_standard

#' * REOMVE THEM
#' 
filtered_prot_dt <- pdt[! PROTEOME_ID %in% c(
    protein_ids_reused, protein_ids_bad_ex2rna, protein_ids_non_standard
    )]

#' Samples still included
unique(filtered_prot_dt$PROTEOME_ID)

#' * Recompute NA frequencies
#' 
compute_na_frequency(filtered_prot_dt, 
    column_id = "GENE_NAME"
)
compute_na_frequency(filtered_prot_dt, 
    column_id = "PROTEOME_ID"
)


#' 
#' ## Save filtered tidy table
#' 
write_tsv(filtered_prot_dt, file=file_tidy_pichler_100min)
stopifnot(file.exists(file_tidy_pichler_100min))


#' 
#' # Pichler QExactive 60min
#' 


