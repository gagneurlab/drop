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
head(grep("uniprot", ignore.case = TRUE, value = TRUE, listAttributes(mart)[,"description"]))
listAttributes(mart)[grep("uniprot", ignore.case = TRUE, listAttributes(mart)[,"description"]),]

# read data 
library(readxl)

norm_prot <- read_excel(file.path(DIR_RAW_PROT, "20180209_robert_proteomics/proteingroups_normalized.xlsx")) %>% as.data.table
setnames(norm_prot, "X__1", "Protein_id")
mat_norm_prot <- as.matrix(t(norm_prot[,2:ncol(norm_prot)]))
colnames(mat_norm_prot) <- norm_prot[,Protein_id]
# TODO check how to filter the samples / genes
table(mat_norm_prot == 0)
mat_norm_prot <- t(mat_norm_prot)
mat_norm_prot[1:10, 1:10]

# look at raw data and correlation

library(LSD)
heatpairs(log10(mat_norm_prot[,c(1:4,18)]))

library(gplots)
heatmap.2(cor(log10(mat_norm_prot + 1)), trace = "none")

# filter prot data
# 
# use only samples where less than 10% of the samples have NA or 0 
sampleZeroFreq <- apply(mat_norm_prot,2,function(x){ sum(x == 0)/length(x) })
table(round(sampleZeroFreq, 3))
mat_norm_prot <- mat_norm_prot[, sampleZeroFreq < 0.1]
dim(mat_norm_prot)
heatmap.2(cor(log10(mat_norm_prot + 1)), trace = "none")

# use only genes where less than of the proteins have 10% NA or 0
protZeroFreq <- apply(mat_norm_prot,1,function(x){ sum(x == 0)/length(x) })
table(round(protZeroFreq, 3))
mat_norm_prot <- mat_norm_prot[protZeroFreq < 0.1,]
dim(mat_norm_prot)
heatmap.2(cor(log10(mat_norm_prot + 1)), trace = "none")

sum(mat_norm_prot < 100)
hist(log10(mat_norm_prot+1), breaks = 100)

# run full APA (aberrant protein abundancies) analysis
#mat_norm_prot <- mat_norm_prot + 1
res <- wrapper_aberrant_protein_expr_simple(prot_intensity=mat_norm_prot, coln_sample_id = "SAMPLEID")

# add HGNC symbol
hgnc_mapping <- data.table(getBM(attributes=c("hgnc_symbol", "uniprotswissprot"),
      filters="uniprotswissprot",
      values=unique(res[,GENE_NAME]),
      mart=mart
))
setnames(res, "GENE_NAME", "UNIPROT_ID")
res <- merge(res, hgnc_mapping, by.x="UNIPROT_ID", by.y="uniprotswissprot", all=TRUE)
setnames(res, "hgnc_symbol", "GENE_NAME")

# check multi mappings
# uniprot id wise
hgnc_mapping[duplicated(uniprotswissprot) | duplicated(uniprotswissprot, fromLast=TRUE)]
# hgnc symbol wise
hgnc_mapping[duplicated(hgnc_symbol) | duplicated(hgnc_symbol, fromLast=TRUE)]


# good results
res[PROT_PADJ < 0.1]

# have fun
hist(res[, PROT_PVALUE])

