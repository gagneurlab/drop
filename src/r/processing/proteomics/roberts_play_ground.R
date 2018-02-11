#' 
#' Roberts Test script
#' 
source("./src/r/config.R")

#' Global variables to use 
#'   * read only [raw]
#'   * processed data to save [proc]
DIR_RAW_PROT  <- "/s/project/mitoMultiOmics/raw_data/proteome"
DIR_PROC_PROT <- "/s/project/mitoMultiOmics/processed_proteomics"

# mapping IDs between gene protein and hgnc_symbol
suppressPackageStartupMessages({
    library(biomaRt)
    library(magrittr)
    library(limma)
    library(DESeq2)
})

?getBM()
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
getBM(attributes = c("affy_hg_u95av2", "hgnc_symbol", "chromosome_name", "band"),
      filters    = "affy_hg_u95av2",
      values     = c("1939_at","1503_at","1454_at"), 
      mart       = mart)

# check out which attribute you need
head(grep("prot", ignore.case = TRUE, value = TRUE, listAttributes(mart)[,"description"]))


# read data 
library(readxl)

norm_prot <- read_excel(file.path(DIR_RAW_PROT, "20180209_robert_proteomics/proteingroups_normalized.xlsx")) %>% as.data.table
setnames(norm_prot, "X__1", "Protein_id")
mat_norm_prot <- as.matrix(t(norm_prot[,2:ncol(norm_prot)]))
colnames(mat_norm_prot) <- norm_prot[,Protein_id]
# TODO check how to filter the samples / genes
table(mat_norm_prot == 0)
mat_norm_prot <- t(mat_norm_prot + 1)
mat_norm_prot[1:10, 1:10]

# look at raw data and correlation

library(LSD)
heatpairs(log10(mat_norm_prot[,c(1:4,18)]))

library(gplots)
heatmap.2(cor(log10(mat_norm_prot)), trace = "none")
heatmap.2(cor(log10(mat_norm_prot[,-18])), trace = "none")

# run full APA (aberrant protein abundancies) analysis
res <- wrapper_aberrant_protein_expr_simple(prot_intensity=mat_norm_prot, coln_sample_id = "SAMPLEID")

# good results
res[PROT_PADJ < 0.1 & SAMPLEID != "P33281 failed"]

# have fun
hist(res[SAMPLEID != "P33281 failed", PROT_PVALUE])

