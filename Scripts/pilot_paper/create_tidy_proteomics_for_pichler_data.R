#'---
#' title: Tidy proteomics for each measurement batch
#' author: Daniel Bader
#' wb:
#'   input: 
#'     - reused_samples: "/s/project/mitoMultiOmics/raw_data//proteome/protein_ids_measured_2nd_time_in_pichler100min.txt"
#'     - bad_exome2rna: "/s/project/mitoMultiOmics/raw_data//proteome/protein_ids_with_bad_exome_to_rna_identity.txt"
#'     - pichler_100min: "/s/project/mitoMultiOmics/raw_data//proteome/20160202_pichler_proteome//pichler100min_combined_corrected_proteinGroups.txt"
#'     - pichler_60min: "/s/project/mitoMultiOmics/raw_data//proteome/20160202_pichler_proteome//pichler60min_combined_corrected_proteinGroups.txt"
#'     - kuester_tmt201706: "/s/project/mitoMultiOmics//raw_data//proteome/20170614_kopajtich_kuester_proteome/m3_lfq_tmt_proteinGroups.txt"
#'   output: 
#'     - raw_pichler_100min: "/s/project/genetic_diagnosis/processed_data/proteome_pichler_100min_unfiltered_all_samples.tsv"
#'     - pichler_100min: "/s/project/genetic_diagnosis/processed_data/proteome_pichler_100min.tsv"
#'     - kuester_tmt201706: "/s/project/genetic_diagnosis/processed_data/proteome_kuester_tmt201706.tsv"
#' output: 
#'   html_document:
#'     toc_float: yes
#'     toc: yes
#'---
#'

#+ echo=F
source("src/r/config.R")

#+ input
file_pichler_100min <- snakemake@input[['pichler_100min']]
file_pichler_60min <- snakemake@input[['pichler_60min']]
# file_pichler_100min  <- "/s/project/mitoMultiOmics/raw_data//proteome/20160202_pichler_proteome//pichler100min_combined_corrected_proteinGroups.txt"
file_kuester_tmt201706 <- "/s/project/mitoMultiOmics//raw_data//proteome/20170614_kopajtich_kuester_proteome/m3_lfq_tmt_proteinGroups.txt"

#+ input data
protein_ids_reused <- readLines(snakemake@input[['reused_samples']])
protein_ids_bad_ex2rna <- readLines(snakemake@input[['bad_exome2rna']])


#' 
#' # Pichler QExactiveHF 100min
#' 

#' 
#' ## Read and process data
#' 
pdt <- wrapper_proteinGroupsTxt_to_tidy_table(
    file_pichler_100min, 
    intensity_column_pattern = "LFQ.intensity."
)
pdt[,ms_method:='Pichler_QExactiveHF_100min']
print(pdt)


#' 
#' ## Save tidy table
#' 
write_tsv(pdt, file= snakemake@output[['raw_pichler_100min']])

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

#' ## Save filtered tidy table
#' 
write_tsv(filtered_prot_dt, file=snakemake@output[['pichler_100min']])



#' 
#' ## Basic statistics
#' 
#' Samples measured
unique(filtered_prot_dt[, NA_FREQ_BY_PROTEOME_ID, by=PROTEOME_ID])

#' Number proteins detected
filtered_prot_dt[PROTEOME_ID=='33254',.N]

#' Avg. NA freq by sample
mean(unique(filtered_prot_dt[, NA_FREQ_BY_PROTEOME_ID, by=PROTEOME_ID])$NA_FREQ_BY_PROTEOME_ID)



#' 
#' # Pichler QExactive 60min
#' 


pdt_pi60 <- wrapper_proteinGroupsTxt_to_tidy_table(
    file_pichler_60min, 
    intensity_column_pattern = "LFQ.intensity."
)
print(pdt_pi60)
pdt_pi60[,ms_method:='Pichler_QExactive_60min']


#' Samples measured
unique(pdt_pi60[, NA_FREQ_BY_PROTEOME_ID, by=PROTEOME_ID])

#' Number proteins detected
pdt_pi60[PROTEOME_ID=='NHDF1',.N]

#' Avg. NA freq by sample
mean(unique(pdt_pi60[, NA_FREQ_BY_PROTEOME_ID, by=PROTEOME_ID])$NA_FREQ_BY_PROTEOME_ID)

#' 
#' ## Summarize technical replicates
#' 

#' ## Save filtered tidy table
#' 


#' 
#' # Kuester TMT labeling from June 2017
#' 

pdt_ku201706 <- wrapper_proteinGroupsTxt_to_tidy_table(
    file_kuester_tmt201706,
    intensity_column_pattern = "Reporter.intensity.corrected."
)
pdt_ku201706[,ms_method:='Kuester_TMT_201706']

#' remove "TMT_10plex_Test_Trinity" samples
pdt_ku201706 <- pdt_ku201706[!grepl("TMT_10plex_Test_Trinity",PROTEOME_ID)]

#' add sample IDs from separate table
file_tmt_sample_map <- file.path(
    'resources', '201706_kuester_tmt_sample_mapping.tsv'
)
tmt_map <- fread(file_tmt_sample_map)
for(i in 0:9){
    pdt_ku201706[PROTEOME_ID==i, PROTEOME_ID:=tmt_map[ID==i, PROTEOME_ID]]
}

print(pdt_ku201706)

#' ## Save filtered tidy table
#' 
write_tsv(pdt_ku201706, file=snakemake@output[['kuester_tmt201706']])


#' 
#' ## Basic statistics
#' 
#' Samples measured
unique(pdt_ku201706[, NA_FREQ_BY_PROTEOME_ID, by=PROTEOME_ID])

#' Number proteins detected
pdt_ku201706[PROTEOME_ID=='nhdf.p9',.N]

#' Avg. NA freq by sample
mean(unique(pdt_ku201706[, NA_FREQ_BY_PROTEOME_ID, by=PROTEOME_ID])$NA_FREQ_BY_PROTEOME_ID)


#'
#' # Merge data sets
#'
print(snakemake@config['rawdir_proteome_kuester_201706'])

