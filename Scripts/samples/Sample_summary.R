#'---
#' title: Samples' summary
#' author: vyepez
#' wb:
#'  input:
#'   - sample_anno: '/s/project/mitoMultiOmics/raw_data/sample_info/SAMPLE_ANNOTATION_PROKISCH.tsv'
#'  output:
#' output:
#'   html_document:
#'    code_folding: show
#'    code_download: TRUE
#'---

#+ echo=F
saveRDS(snakemake, "tmp/sample_summary.snakemake")
# snakemake <- readRDS("tmp/sample_summary.snakemake")

suppressPackageStartupMessages({
    library(ggplot2)
    library(ggthemes)
    library(data.table)
    library(magrittr)
    library(dplyr)
})

col_unsolved = "gray75"
col_solved = "mediumseagreen"

sa <- fread(snakemake@input$sample_anno)

#' Number of patient fibroblasts in glucose with both RNA and WES
sa_pat <- sa[!is.na(RNA_ID) & !is.na(EXOME_ID) & TISSUE == 'FIBROBLAST' & GROWTH_MEDIUM == 'GLU' & DISEASE != 'HEALTHY' & is.na(TRANSDUCED_GENE) & FIBROBLAST_ID != 'NHDF'] 
uniqueN(sa_pat[, PATIENT_ID])

DT::datatable(sa_pat[, .(PATIENT_ID, FIBROBLAST_ID, EXOME_ID, RNA_ID, KNOWN_MUTATION, BATCH, DISEASE, COMMENT)])

#' ## Samples without patient it
DT::datatable(sa_pat[is.na(PATIENT_ID)])

#' ## Duplicated samples
dup_patients <- sa_pat[duplicated(PATIENT_ID), unique(PATIENT_ID)] 
DT::datatable(sa_pat[!is.na(PATIENT_ID)][PATIENT_ID %in% dup_patients])

#' ## Read and clean global sample annotation
sa <- sa[!is.na(BATCH)]
sa[, solved := !is.na(KNOWN_MUTATION)]
sa[, sample_type := 'patient']
sa[DISEASE == 'HEALTHY' | FIBROBLAST_ID == 'NHDF', sample_type := 'control']
sa[!is.na(TRANSDUCED_GENE) | TISSUE != 'FIBROBLAST' | GROWTH_MEDIUM != 'GLU', sample_type := 'other']
sa[, IS_RNA_SEQ_STRANDED := as.logical(IS_RNA_SEQ_STRANDED)]

ggplot(sa, aes(BATCH)) + geom_bar(aes(y = ..count.., fill = sample_type)) + 
    geom_text(aes(label = ..count..), stat = 'count' , vjust = -.5) + 
    theme_bw(base_size = 14) + scale_fill_canva(palette="Subdued and proffesional")

ggplot(sa, aes(BATCH)) + geom_bar(aes(y = ..count.., fill = RNA_PERSON)) + 
    theme_bw(base_size = 14) + scale_fill_ptol()

#' Total number of samples coming from patients (strand and non strand specific)
sa[sample_type == 'patient', .N]
sa[sample_type == 'patient', table(IS_RNA_SEQ_STRANDED)]   

ggplot(sa[sample_type == 'patient'], aes(BATCH)) + geom_bar(aes(y = ..count.., fill = solved)) + 
    geom_text(aes(label = ..count..), stat = 'count' , vjust = -.5) + 
    theme_bw(base_size = 14) + scale_fill_manual(values = c(col_unsolved, col_solved)) + ggtitle('Patients only')

ggplot(sa[sample_type == 'patient'], aes(IS_RNA_SEQ_STRANDED)) + geom_bar(aes(y = ..count.., fill = solved)) + 
    geom_text(aes(label = ..count..), stat = 'count' , vjust = -.5) + 
    theme_bw(base_size = 14) + scale_fill_manual(values = c(col_unsolved, col_solved)) + ggtitle('Patients only')

#' ## Download sample annotation with patients only
#' ### Add Sarah's diagnose
exo_db <- fread("/s/project/mitoMultiOmics/exomes1000/processed_data/exomes_clean.tsv")
exo_db = exo_db[WES_in_db == 'Y']
sa2 <- left_join(sa, exo_db[, .(RNA_ID, gene_name)]) %>% as.data.table
sa2[gene_name == '.', gene_name := NA]

#' ### Diagnose mismatches
DT::datatable(sa2[KNOWN_MUTATION != gene_name])
setnames(sa2, 'gene_name', 'sarah_diagnose')

#' ### Download table
#+ echo=F
write.table(sa2[sample_type == 'patient'][order(BATCH)], "/s/public_webshare/project/genetic_diagnosis/results/sample_anno_rna_patients.txt", sep = "\t", quote = F, row.names = F)
#' [Download rna sample annotation](https://i12g-gagneurweb.informatik.tu-muenchen.de/project/genetic_diagnosis/results/sample_anno_rna_patients.txt)
#'

#' ## Questions to answer
#' ### 1. Sample 62336R appears both in Batch 2 and Batch 4, both are repetition from Laura; which should be included?
#' ### 2. What is the real diagnosis of EXOME 65125?


#' ## Check exomes in our Gagneurlab server and save them
DIR_raw <- "/s/project/mitoMultiOmics/raw_data/helmholtz/"
list_files <- sapply(list.files(DIR_raw), function(f) list.files(file.path(DIR_raw, f)))
exomes_our <- names(list_files["exomicout" %in% list_files])
exomes_not_our_db <- sort(setdiff(sa$EXOME_ID, exomes_our))
length(exomes_not_our_db)
write.table(exomes_not_our_db, "resources/missing_exomes.txt", sep = "\n", row.names = F, col.names = F, quote = F)
