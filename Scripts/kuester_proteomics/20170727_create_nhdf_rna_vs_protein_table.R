#'---
#' title: Compare RNA vs Protein for NHDF
#' author: Daniel Bader
#' wb:
#'   input: [ 
#'   "/s/project/mitoMultiOmics//raw_data//proteome/20170614_kopajtich_kuester_proteome/m3_lfq_tmt_proteinGroups.txt",
#'   "/s/project/mitoMultiOmics//raw_data//proteome/20170614_kopajtich_kuester_proteome/m4_lfq_id_trinity_proteinGroups.txt",
#'   "/s/project/mitoMultiOmics/counttable_galignment/rna/raw_counttable_rna_no_strand.tsv",
#'   "resources/201706_kuester_tmt_sample_mapping.tsv"
#'   ]
#'   output: "/s/project/genetic_diagnosis/processed_data/nhdf_rna_proteome_tmt_trinity.tsv"
#' output: 
#'   html_document:
#'     toc_float: yes
#'     toc: yes
#'---
#'



#+ echo=F
source("src/r/config.R")
#+ 
opts_chunk$set(message=T)

#+ in out data
protdir <- file.path(RAWDIR, "proteome", "20170614_kopajtich_kuester_proteome")
files_kuester <- list.files(protdir, pattern = "^m.*txt$", full.names = T)
file_rna_raw_counts <- file.path(
    '/s/project/mitoMultiOmics/counttable_galignment/rna/raw_counttable_rna_no_strand.tsv'
)
file_tmt_sample_map <- file.path(
    'resources', '201706_kuester_tmt_sample_mapping.tsv'
)

# out
file_nhdf_out <- file.path(PROC_DATA, "nhdf_rna_proteome_tmt_trinity.tsv")

#' 
#' # Get iBAQ proteome expression
#' 

#' ## TMT   
pdt_tmt <- wrapper_proteinGroupsTxt_to_tidy_table(
    files_kuester[3],
    intensity_column_pattern = "Reporter.intensity.corrected.",
    column_intensity = "LFQ_INTENSITY_TMT"
)

#' remove "TMT_10plex_Test_Trinity" samples
pdt_tmt <- pdt_tmt[!grepl("TMT_10plex_Test_Trinity",PROTEOME_ID)]

#' add sample IDs from separate table
tmt_map <- fread(file_tmt_sample_map)
for(i in 0:9){
    pdt_tmt[PROTEOME_ID==i, PROTEOME_ID:=tmt_map[ID==i, PROTEOME_ID]]
}



#' ## Trinity
pdt_trinity <- wrapper_proteinGroupsTxt_to_tidy_table(
    files_kuester[4],
    intensity_column_pattern = "Intensity.",
    column_intensity = "LFQ_INTENSITY_TRINITY"
)
pdt_trinity[, PROTEOME_ID:="nhdf.p9"]


pdt_trinity_ibaq <- wrapper_proteinGroupsTxt_to_tidy_table(
    files_kuester[4],
    intensity_column_pattern = "iBAQ\\.",
    column_intensity = "IBAQ_TRINITY"
)
pdt_trinity_ibaq[, PROTEOME_ID:="nhdf.p9"]


#' 
#' ## merge
#' 

#' * reformat tables
nhdf_tmt <- dcast(
    pdt_tmt[grepl("nhdf", PROTEOME_ID)], 
    GENE_NAME ~ PROTEOME_ID, 
    value.var = 'LFQ_INTENSITY_TMT'
    )
setnames(nhdf_tmt, gsub('nhdf', 'lfq_tmt_nhdf', names(nhdf_tmt)))


nhdf_trinity_lfq <- dcast(pdt_trinity, 
    GENE_NAME ~ PROTEOME_ID, 
    value.var = 'LFQ_INTENSITY_TRINITY'
)
setnames(nhdf_trinity_lfq, gsub('nhdf', 'lfq_trinity_nhdf', names(nhdf_trinity_lfq)))


nhdf_trinity_ibaq <- dcast(pdt_trinity_ibaq, 
    GENE_NAME ~ PROTEOME_ID, 
    value.var = 'IBAQ_TRINITY'
)
setnames(nhdf_trinity_ibaq, gsub('nhdf', 'ibaq_trinity_nhdf', names(nhdf_trinity_ibaq)))

#' * combine
#' 
nhdf_proteome <- merge(
    merge(nhdf_tmt, nhdf_trinity_lfq, all=TRUE), 
    nhdf_trinity_ibaq, all=TRUE
)


#' 
#' # Get RNA NHDF raw counts
#' 

#' * get NHDF RNA ids
#' 
rna_ids_nhdf <- get_rna_for_fibro(
    'NHDF', 
    SAMPLE_ANNOTATION[GROWTH_MEDIUM=='GLU' & is.na(TRANSDUCED_GENE) & TISSUE=='FIBROBLAST']
    )
rna_ids_nhdf
rna_ids_nhdf = rna_ids_nhdf[1:8]

#' * read and subset raw counts
#' 
rna_res <- as.data.table(read.delim(file_rna_raw_counts, check.names = F), keep.rownames=T)
setnames(rna_res, 'rn', 'GENE_NAME')
rna_nhdf <- rna_res[, c('GENE_NAME',rna_ids_nhdf), with=F]
setnames(rna_nhdf, gsub('^[^G]', 'rna_gene_counts_nhdf_', names(rna_nhdf)))


#'
#' # merge rna protein
#' 

nhdf_combo <- merge(nhdf_proteome, rna_nhdf, all=TRUE)
write_tsv(nhdf_combo, file = file_nhdf_out)



