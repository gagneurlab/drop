#'---
#' title: Proteomics Samples Overview
#' author: vyepez
#' wb:
#'  input: 
#'   - protein_sa: '/s/project/mitoMultiOmics/raw_data/proteome/20180927_robert_proteomics/raw_data/P014_10_template_for_row_wise_normalization_BBM_overview.xlsx'
#'   - protein_mat: '/s/project/mitoMultiOmics/raw_data/proteome/20180927_robert_proteomics/TMT_1_14_row_col_norm_all_wo_dup_wo_contr.xlsx'
#' output: 
#'   html_document:
#'    code_folding: show
#'    code_download: TRUE
#'---

#+ echo=F
saveRDS(snakemake, "tmp/prot_overview.snakemake")
# snakemake <- readRDS("tmp/prot_overview.snakemake")
suppressPackageStartupMessages({
    library(readxl)
    library(data.table)
    library(magrittr)
    library(biomaRt)
    library(ggplot2)
    library(ggthemes)
    library(dplyr)
})

#' ## Read and clean both tables
#' ### Read sample annotation
sa_prot <- read_xlsx(snakemake@input$protein_sa) %>% as.data.table
setnames(sa_prot, "ID sample", "PROTEOME_ID")
sa_prot[, Use := TRUE]
sa_prot[grep("errorenous sample|control sample|technical replicates", Comment), Use := FALSE]

sa_prot[Use == T, s_run := 1:.N, by = Run]
sa_prot[, Batch := substr(gsub("mix", "", Run), 1, 1)]
sa_prot[Use == T, s_batch := 1:.N, by = Batch]
sa_prot[, Subgroup := substr(gsub("mix", "", Run), 2, 2)]

#' Valid samples in sample annotation (w/o controls or replicates)
sum(sa_prot$Use)


#' ### Read protein intensities matrix
mat_all <- read_xlsx(snakemake@input$protein_mat) %>% as.data.table
dim(mat_all)

# remove outlier samples
mat_all <- mat_all[, -c("P33281", "P103213")]
dim(mat_all)

# check if proteins is 0 in all samples (yes if TRUE)
any(rowSums(mat_all[, -"ProteinID"]) == 0)
sum(rowSums(mat_all[, -"ProteinID"]) == 0)

# remove rows with all zeroes
mat_all <- mat_all[! rowSums(mat_all[, -"ProteinID"]) == 0]

protein_mat <- as.matrix(mat_all[, -"ProteinID"])
row.names(protein_mat) <- mat_all$ProteinID
samples_take <- sa_prot[s_batch == 1, PROTEOME_ID]
upset_mat <- (protein_mat[, samples_take] > 0) * 1
upset_mat[is.na(upset_mat)] <- 0 
colnames(upset_mat) <- sa_prot[s_batch == 1, Batch]

library(UpSetR)
upset_dt <- data.table(Protein_ID = mat_all$ProteinID, upset_mat)
upset(upset_dt)

#' Proteins detected in at least 1 sample
nrow(mat_all)

samples_take <- sa_prot[s_run == 1, PROTEOME_ID]
detected_prot <- data.table(PROTEOME_ID = colnames(protein_mat[, samples_take]), Proteins_detected = colSums(protein_mat[,samples_take]>0))
detected_prot <- left_join(detected_prot, sa_prot[,.(PROTEOME_ID, Batch, Subgroup)]) %>% as.data.table
library(ggbeeswarm)
#+ 
ggplot(detected_prot, aes(Batch, Proteins_detected, fill = Subgroup)) + 
    geom_bar(stat = 'identity', position = 'dodge') + 
    geom_text(aes(label = Proteins_detected), vjust = -.3, stat = 'identity', position = position_dodge(width=.9)) +
    theme_bw() + scale_fill_canva() + coord_cartesian(ylim=c(6e3,8e3))

#' ## Map genes to proteins

# create BiomaRt object for mapping
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

# add HGNC symbol
hgnc_mapping <- data.table(getBM(attributes=c("hgnc_symbol", "uniprotswissprot"),
                                 filters="uniprotswissprot",
                                 values=row.names(protein_mat),
                                 mart=mart
))
# setnames(res, "GENE_NAME", "UNIPROT_ID")
# res <- merge(res, hgnc_mapping, by.x="UNIPROT_ID", by.y="uniprotswissprot", all=TRUE)
# setnames(res, "hgnc_symbol", "GENE_NAME")

# check multi mappings
# uniprot id wise
hgnc_mapping[duplicated(uniprotswissprot) | duplicated(uniprotswissprot, fromLast=TRUE)]
# hgnc symbol wise
hgnc_mapping[duplicated(hgnc_symbol) | duplicated(hgnc_symbol, fromLast=TRUE)]

