#'---
#' title: Aberrant protein expression using LIMMA
#' author: Daniel Bader
#' wb:
#'   input: [
#'     proteome_merged: "/s/project/genetic_diagnosis/processed_data/proteome_intensities_all_merged.tsv"
#'   ]
#'   output: [
#'     proteome_aberexp: "/s/project/genetic_diagnosis/processed_results/proteome_aberrant_expression.tsv"
#'   ]
#' output: 
#'   html_document:
#'     toc_float: yes
#'     toc: yes
#'---
#'

#+ echo=F
source('.wBuild/wBuildParser.R')
parseWBHeader("Scripts/proteomics/create_aberrant_protein_expression.R")

source("src/r/config.R")
file_proteome_merged <- snakemake@input[['proteome_merged']]
file_aber_prot_exp <- snakemake@output[['proteome_aberexp']]


#' Input table tidy raw proteome LFQ intensities
proteome_merged <- fread(file_proteome_merged)
head(proteome_merged)

names(proteome_merged) = toupper(names(proteome_merged))

#' Available proteomics measurements:
all_ms_methods <- unique(proteome_merged$MS_METHOD)
all_ms_methods


#'
#' # Aberrant expression
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

# result object
prot_aberexp <- data.table()

#' Perform test separately for each measurement
#' 
for(tmp_ms_method in all_ms_methods){
    message(tmp_ms_method)
    # select one measurement
    prot_intensity_dt <- proteome_merged[MS_METHOD==tmp_ms_method]
    # get matrix representation for limma
    prot_intensity_mat <- convert_tidy_table_to_numeric_matrix(
        prot_intensity_dt, 'GENE_NAME', 'PROTEOME_ID', 'LFQ_INTENSITY'
    )
    
    # differential expression with limma
    single_prot_aberexp <- wrapper_aberrant_protein_expr_simple(prot_intensity_mat)
    # compute Z-score with fold change from limma 
    # and sd of normalized log intensities
    single_prot_aberexp <- compute_protein_zscore(single_prot_aberexp)
    
    # merge with raw intensities
    setnames(prot_intensity_dt, 'PROTEOME_ID', 'FIBROBLAST_ID')
    single_prot_aberexp <- merge(
        prot_intensity_dt, 
        single_prot_aberexp, 
        by=c('FIBROBLAST_ID', 'GENE_NAME'), 
        all=T
    )
    
    # add to result object
    prot_aberexp <- rbind(prot_aberexp, single_prot_aberexp)
}

#' ## Examine output
head(prot_aberexp)

#' Number of samples by Method
unique(prot_aberexp[, FIBROBLAST_ID, by=MS_METHOD])[, .N, by=MS_METHOD]
#' Number of proteins by Method
unique(prot_aberexp[, .N, by=c('FIBROBLAST_ID', 'MS_METHOD')][,.(N, MS_METHOD)])


#' SAVE
file_aber_prot_exp
write_tsv(prot_aberexp, file = file_aber_prot_exp)


#' ## Visualize
res_prot <- fread("/s/project/genetic_diagnosis/processed_results/proteome_aberrant_expression.tsv")
dim(res_prot)
#' number of samples
res_prot[, uniqueN(FIBROBLAST_ID)]
res_prot[, FC := round(2^PROT_LOG2FC, 3)]
res_prot[, c("NA_FREQ_BY_GENE_NAME", "NA_FREQ_BY_PROTEOME_ID", "PROT_LOGODDS", "SD_NORM_LOG2_LFQ") := NULL]
# res_prot[FIBROBLAST_ID == 62335 & GENE_NAME == 'NFU1']
DT::datatable(res_prot[PROT_PADJ < .1], filter = 'top')

