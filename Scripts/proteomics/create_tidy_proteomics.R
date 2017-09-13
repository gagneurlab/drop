#'---
#' title: Tidy proteomics for each measurement batch
#' author: Daniel Bader
#' wb:
#'   input: 
#'     - reused_samples: "`sm config['rawdir_proteome'] + 'protein_ids_measured_2nd_time_in_pichler100min.txt'`"
#'     - bad_exome2rna: "`sm config['rawdir_proteome'] + 'protein_ids_with_bad_exome_to_rna_identity.txt'`"
#'     - pichler_100min: "`sm config['rawdir_proteome'] + '20160202_pichler_proteome//pichler100min_combined_corrected_proteinGroups.txt'`"
#'     - pichler_60min: "`sm config['rawdir_proteome'] + '20160202_pichler_proteome//pichler60min_combined_corrected_proteinGroups.txt'`"
#'     - kuester_tmt201706: "`sm config['rawdir_proteome'] + '20170614_kopajtich_kuester_proteome/m3_lfq_tmt_proteinGroups.txt'`"
#'   output: 
#'     - raw_pichler_100min: "`sm config['PROC_DATA'] + 'proteome_pichler_100min_unfiltered_all_samples.tsv'`"
#'     - pichler_100min: "`sm config['PROC_DATA'] + 'proteome_pichler_100min.tsv'`"
#'     - raw_pichler_60min: "`sm config['PROC_DATA'] + 'proteome_pichler_60min_unfiltered_all_samples.tsv'`"
#'     - pichler_60min: "`sm config['PROC_DATA'] + 'proteome_pichler_60min.tsv'`"
#'     - raw_kuester_tmt201706: "`sm config['PROC_DATA'] + 'proteome_kuester_tmt201706_unfiltered_all_samples.tsv'`"
#'     - kuester_tmt201706: "`sm config['PROC_DATA'] + 'proteome_kuester_tmt201706.tsv'`"
#'     - tidy_proteomics: "`sm config['PROC_DATA'] + 'proteome_intensities_all_merged.tsv'`"
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
file_kuester_tmt201706 <- snakemake@input[['kuester_tmt201706']]
file_kuester_tmt_sample_map <- file.path(
    'resources', '201706_kuester_tmt_sample_mapping.tsv'
)
#+ input data
protein_ids_reused <- readLines(snakemake@input[['reused_samples']])
protein_ids_bad_ex2rna <- readLines(snakemake@input[['bad_exome2rna']])


#' 
#' # Pichler QExactiveHF 100min
#' 

#' 
#' ## Read and process data
#' 
raw_prot_pi100 <- wrapper_proteinGroupsTxt_to_tidy_table(
    file_pichler_100min, 
    intensity_column_pattern = "LFQ.intensity."
)
raw_prot_pi100[,MS_METHOD:='Pichler_QExactiveHF_100min']
head(raw_prot_pi100)


#' 
#' ## Save tidy table
#' 
write_tsv(raw_prot_pi100, file= snakemake@output[['raw_pichler_100min']])

#' Samples measured
unique(raw_prot_pi100$PROTEOME_ID)


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
    grep("2w", unique(raw_prot_pi100$PROTEOME_ID), value=T),
    grep("\\.T\\.", unique(raw_prot_pi100$PROTEOME_ID), value=T)
)
protein_ids_non_standard

#' * REOMVE THEM
#' 
filtered_prot_pi100 <- raw_prot_pi100[! PROTEOME_ID %in% c(
    protein_ids_reused, protein_ids_bad_ex2rna, protein_ids_non_standard
    )]

#' Samples still included
unique(filtered_prot_pi100$PROTEOME_ID)

#' * Recompute NA frequencies
#' 
compute_na_frequency(filtered_prot_pi100, 
    column_id = "GENE_NAME"
)
compute_na_frequency(filtered_prot_pi100, 
    column_id = "PROTEOME_ID"
)

#' ## Save filtered tidy table
#' 
write_tsv(filtered_prot_pi100, file=snakemake@output[['pichler_100min']])




#' 
#' # Pichler QExactive 60min
#' 
#file_pichler_60min <- "/s/project/mitoMultiOmics/raw_data//proteome/20160202_pichler_proteome//pichler60min_combined_corrected_proteinGroups.txt"

raw_pdt_pi60 <- wrapper_proteinGroupsTxt_to_tidy_table(
    file_pichler_60min, 
    intensity_column_pattern = "LFQ.intensity."
)
raw_pdt_pi60[,MS_METHOD:='Pichler_QExactive_60min']
head(raw_pdt_pi60)

write_tsv(raw_pdt_pi60, file= snakemake@output[['raw_pichler_60min']])

#' 
#' ## Summarize technical replicates
#' 

#' Build mean over technical replicates ignoring missing values. 
#' Which is the same as imputing them with the mean, 
#' i.e. the other value from one of the 2 replicates.

raw_pdt_pi60[, sample:=gsub('*_..','', PROTEOME_ID)]
raw_pdt_pi60[grepl('NHDF', PROTEOME_ID), sample:='NHDF']
# reduce by mean
reduced_pi60 <- raw_pdt_pi60[, .
    (LFQ_INTENSITY=mean(LFQ_INTENSITY, na.rm=T)), 
    by=c("sample", "PROTEIN_ID", "GENE_NAME")
]
reduced_pi60[,MS_METHOD:=unique(raw_pdt_pi60$MS_METHOD)]
setnames(reduced_pi60, 'sample', 'PROTEOME_ID')

#' * Recompute NA frequencies
#' 
compute_na_frequency(reduced_pi60, 
    column_id = "GENE_NAME"
)
compute_na_frequency(reduced_pi60, 
    column_id = "PROTEOME_ID"
)


#' ## Save filtered tidy table
#' 
write_tsv(reduced_pi60, file= snakemake@output[['pichler_60min']])




#' 
#' # Kuester TMT labeling from June 2017
#' 

raw_pdt_ku201706 <- wrapper_proteinGroupsTxt_to_tidy_table(
    file_kuester_tmt201706,
    intensity_column_pattern = "Reporter.intensity.corrected."
)
raw_pdt_ku201706[,MS_METHOD:='Kuester_TMT_201706']

#' * remove "TMT_10plex_Test_Trinity" samples, 
#'   since they are duplicated columns.
raw_pdt_ku201706 <- raw_pdt_ku201706[!grepl("TMT_10plex_Test_Trinity",PROTEOME_ID)]

#' * add sample IDs from separate table
tmt_map <- fread(file_kuester_tmt_sample_map)
for(i in 0:9){
    raw_pdt_ku201706[PROTEOME_ID==i, PROTEOME_ID:=tmt_map[ID==i, PROTEOME_ID]]
}

head(raw_pdt_ku201706)
unique(raw_pdt_ku201706$PROTEOME_ID)
write_tsv(raw_pdt_ku201706, file=snakemake@output[['raw_kuester_tmt201706']])



#' 
#' ## Summarize technical replicates
#' 

#' Build mean over technical replicates ignoring missing values. 
#' Which is the same as imputing them with the mean, 
#' i.e. the other value from one of the 2 replicates.

raw_pdt_ku201706[, sample:=gsub('(*)\\..*$', toupper('\\1'), PROTEOME_ID)]
raw_pdt_ku201706[, sample:= toupper(sample)]
unique(raw_pdt_ku201706$sample)
# reduce by mean
reduced_pdt_ku201706 <- raw_pdt_ku201706[, 
    .(LFQ_INTENSITY=mean(LFQ_INTENSITY, na.rm=T)), 
    by=c("sample", "PROTEIN_ID", "GENE_NAME")]
reduced_pdt_ku201706[,MS_METHOD:=unique(raw_pdt_ku201706$MS_METHOD)]
setnames(reduced_pdt_ku201706, 'sample', 'PROTEOME_ID')

#' * Recompute NA frequencies
#' 
compute_na_frequency(reduced_pdt_ku201706, 
    column_id = "GENE_NAME"
)
compute_na_frequency(reduced_pdt_ku201706, 
    column_id = "PROTEOME_ID"
)



#' ## Save filtered tidy table
#' 
write_tsv(reduced_pdt_ku201706, file=snakemake@output[['kuester_tmt201706']])



#'
#' # Merge data sets
#'

proteomics_tidy_dt <- rbind(
    filtered_prot_pi100, reduced_pi60, reduced_pdt_ku201706
)
head(proteomics_tidy_dt)

write_tsv(proteomics_tidy_dt, file=snakemake@output[['tidy_proteomics']])

#' 
#' ## Basic statistics
#' 
#' Samples measured
unique(proteomics_tidy_dt[, PROTEOME_ID, by=MS_METHOD])[,.N, by=MS_METHOD]

#' Number proteins detected
unique(proteomics_tidy_dt[, .N, by=c("PROTEOME_ID","MS_METHOD")][,.(MS_METHOD, N)])

#' Avg. NA freq by sample
proteomics_tidy_dt[, 
    .(mean_NA_freq_by_sample=mean(NA_FREQ_BY_PROTEOME_ID)), 
    by=MS_METHOD
]

