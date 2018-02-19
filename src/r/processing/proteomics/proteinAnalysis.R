
#
# 2018_02_16 TEST
#

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

# read data 
library(readxl)

norm_prot <- read_excel(file.path(DIR_RAW_PROT, "20180209_robert_proteomics/proteingroups_rowwise_normalized_with_zeroblocks.xlsx")) %>% as.data.table
setnames(norm_prot, "X__1", "ProteinID")
mat_norm_prot <- as.matrix(t(norm_prot[,2:ncol(norm_prot)]))
colnames(mat_norm_prot) <- norm_prot[,ProteinID]

# look at raw data and correlation

library(LSD)
heatpairs(log10(mat_norm_prot[,c(1:4,18)]))

library(gplots)
heatmap.2(cor(log10(mat_norm_prot + 1)), trace = "none")as.data.table
setnames(norm_prot_zero, "X__1", "ProteinID")
mat_norm_prot_zero <- as.matrix(t(norm_prot_zero[,2:ncol(norm_prot_zero)]))
colnames(mat_norm_prot_zero) <- norm_prot_zero[,ProteinID]


# look at raw data and correlation

library(LSD)
heatpairs(log10(mat_norm_prot_zero[,c(1:4,18)]))

library(gplots)
heatmap.2(cor(log10(mat_norm_prot_zero + 1)), trace = "none")
