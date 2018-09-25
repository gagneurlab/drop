#' 
#' Roberts Test script
#' 
source("./src/r/config.R")
library(R.utils)
sourceDirectory("../gagneurlab_shared/r/knitr_helper/")

#' Global variables to use 
#'  * read only [raw]
#'  * processed data to save [proc]
DIR_RAW_PROT  <- "/s/project/mitoMultiOmics/raw_data/proteome"
DIR_PROC_PROT <- "/s/project/mitoMultiOmics/processed_proteomics"

# mapping IDs between gene protein and hgnc_symbol
suppressPackageStartupMessages({
    library(biomaRt)
    library(magrittr)
    library(limma)
    library(DESeq2)
    library(tidyr)
    library(readxl)
    library(LSD)
    library(gplots)
    library(knitr)
    library(markdown)
})

# create BiomaRt object for mapping
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

# render_webserver_html("./src/r/processing/proteomics/roberts_play_ground.R",project_webserver_dir = "/s/public_webshare/project/genetic_diagnosis/proteome/")


# read data 
norm_prot_zero <- read_excel(file.path(DIR_RAW_PROT, "20180209_robert_proteomics/proteingroups_rowwise_normalized_with_zeroblocks.xlsx")) %>% as.data.table
mat_norm_prot_zero <- as.matrix(norm_prot_zero[,2:ncol(norm_prot_zero)])
rownames(mat_norm_prot_zero) <- norm_prot_zero[,ProteinID]

# TODO check how to filter the samples / genes
table(mat_norm_prot_zero == 0)
mat_norm_prot_zero[1:10, 1:10]

# look at raw data and correlation
# heatpairs(log10(mat_norm_prot_zero[,c(1:4,18)]))
heatmap.2(cor(log10(mat_norm_prot_zero + 1)), trace = "none")

#' 
#' # Filter prot data
# 
# use only samples with less than 10% of NA or 0 
sampleZeroFreq <- apply(mat_norm_prot_zero,2,function(x){ sum(x == 0)/length(x) })
table(round(sampleZeroFreq, 3))
mat_norm_prot_zero <- mat_norm_prot_zero[, sampleZeroFreq < 0.1]
dim(mat_norm_prot_zero)
heatmap.2(cor(log10(mat_norm_prot_zero + 1)), trace = "none")

# use only genes where less than of the proteins have 10% NA or 0
protZeroFreq <- apply(mat_norm_prot_zero,1,function(x){ sum(x == 0)/length(x) })
table(round(protZeroFreq, 3))
mat_norm_prot_zero <- mat_norm_prot_zero[protZeroFreq < 0.1,]
dim(mat_norm_prot_zero)
heatmap.2(cor(log10(mat_norm_prot_zero + 1)), trace = "none")

sum(mat_norm_prot_zero < 100)
hist(log10(mat_norm_prot_zero+1), breaks = 100)

# run full APA (aberrant protein abundancies) analysis
#mat_norm_prot <- mat_norm_prot + 1
res <- wrapper_aberrant_protein_expr_simple(prot_intensity=mat_norm_prot_zero, coln_sample_id = "SAMPLEID")

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
res[PROT_PADJ < 0.1][order(PROT_PADJ)]

# have fun
hist(res[, PROT_PVALUE])


# calc zscore
l2fc_mat <- t(matrix(res[,PROT_LOG2FC], nrow=ncol(mat_norm_prot_zero)))
zscores <- (l2fc_mat - rowMeans(l2fc_mat))/rowSds(l2fc_mat)
hist(zscores)
table(abs(zscores) >= 3)

res[,L2FC_ZSCORE:=as.vector(zscores)]

res[, signif := PROT_PADJ < .05]
res[, N := sum(signif, na.rm = T), by = SAMPLEID]

res_sub = res[signif == TRUE]
View(res_sub[N < 3])

setorder(res_sub, N)


#' 
#' # Results
#' 
DT::datatable(res[order(PROT_PADJ)])

# write it out to disk
write_tsv(res, file.path(DIR_PROC_PROT, "20180219_myResultTable.tsv"))









#
# TEST
#
TMT3_zero <- read_excel(file.path(DIR_RAW_PROT, "20180209_robert_proteomics/proteingroups_rowwise_normalized_with_zeroblocks_TMT3.xlsx")) %>% as.data.table
mat_TMT3_zero <- as.matrix(TMT3_zero[,2:ncol(TMT3_zero)])
rownames(mat_TMT3_zero) <- TMT3_zero[,ProteinID]

# TODO check how to filter the samples / genes
table(mat_TMT3_zero == 0)
mat_TMT3_zero[1:10, 1:10]

# look at raw data and correlation
# heatpairs(log10(mat_TMT1_zero[,c(1:4,18)]))
heatmap.2(cor(log10(mat_TMT3_zero + 1)), trace = "none")

#' 
#' # Filter prot data
# 
# use only samples with less than 10% of NA or 0 
sampleZeroFreq <- apply(mat_TMT3_zero,2,function(x){ sum(x == 0)/length(x) })
table(round(sampleZeroFreq, 3))
mat_TMT3_zero <- mat_TMT3_zero[, sampleZeroFreq < 0.1]
dim(mat_TMT3_zero)
heatmap.2(cor(log10(mat_TMT3_zero + 1)), trace = "none")

# use only genes where less than of the proteins have 10% NA or 0
protZeroFreq <- apply(mat_TMT3_zero,1,function(x){ sum(x == 0)/length(x) })
table(round(protZeroFreq, 3))
mat_TMT3_zero <- mat_TMT3_zero[protZeroFreq < 0.1,]
dim(mat_TMT3_zero)
heatmap.2(cor(log10(mat_TMT3_zero + 1)), trace = "none")

sum(mat_TMT3_zero < 100)
hist(log10(mat_TMT3_zero+1), breaks = 100)

# run full APA (aberrant protein abundancies) analysis
#mat_TMT1_zero <- mat_TMT1_zero + 1
res <- wrapper_aberrant_protein_expr_simple(prot_intensity=mat_TMT3_zero, coln_sample_id = "SAMPLEID")

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
res[PROT_PADJ < 0.1][order(PROT_PADJ)]

# have fun
hist(res[, PROT_PVALUE])


# calc zscore
l2fc_mat_TMT3 <- t(matrix(res[,PROT_LOG2FC], nrow=ncol(mat_TMT3_zero)))
zscores <- (l2fc_mat_TMT3 - rowMeans(l2fc_mat))/rowSds(l2fc_mat)
hist(zscores, breaks=50)
table(abs(zscores) >= 3)

res[,L2FC_ZSCORE:=as.vector(zscores)]

#' 
#' # Results
#' 
DT::datatable(res[order(PROT_PADJ)])

# write it out to disk
write_tsv(res, file.path(DIR_PROC_PROT, "20180226_myResultTable_TMT3.tsv"))



#
# TEST on old MS2 TMT dataset
#
TMT_MS2 <- read_excel(file.path(DIR_RAW_PROT, "20180209_robert_proteomics/TMT_MS2.xlsx")) %>% as.data.table
mat_TMT_MS2 <- as.matrix(TMT_MS2[,2:ncol(TMT_MS2)])
rownames(mat_TMT_MS2) <- TMT_MS2[,ProteinID]

# TODO check how to filter the samples / genes
table(mat_TMT_MS2 == 0)
mat_TMT_MS2[1:10, 1:9]

# look at raw data and correlation
# heatpairs(log10(mat_TMT_MS2[,c(1:4,18)]))
heatmap.2(cor(log10(mat_TMT_MS2 + 1)), trace = "none")

#' 
#' # Filter prot data
# 
# use only samples with less than 10% of NA or 0 
sampleZeroFreq <- apply(mat_TMT_MS2,2,function(x){ sum(x == 0)/length(x) })
table(round(sampleZeroFreq, 3))
mat_TMT_MS2 <- mat_TMT_MS2[, sampleZeroFreq < 0.1]
dim(mat_TMT_MS2)
heatmap.2(cor(log10(mat_TMT_MS2 + 1)), trace = "none")

# use only genes where less than of the proteins have 10% NA or 0
protZeroFreq <- apply(mat_TMT_MS2,1,function(x){ sum(x == 0)/length(x) })
table(round(protZeroFreq, 3))
mat_TMT_MS2 <- mat_TMT_MS2[protZeroFreq < 0.1,]
dim(mat_TMT_MS2)
heatmap.2(cor(log10(mat_TMT_MS2 + 1)), trace = "none")

sum(mat_TMT_MS2 < 100)
hist(log10(mat_TMT_MS2+1), breaks = 100)

# run full APA (aberrant protein abundancies) analysis
#mat_TMT1_zero <- mat_TMT1_zero + 1
res <- wrapper_aberrant_protein_expr_simple(prot_intensity=mat_TMT_MS2, coln_sample_id = "SAMPLEID")

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
res[PROT_PADJ < 0.1][order(PROT_PADJ)]

# have fun
hist(res[, PROT_PVALUE])


# calc zscore
l2fc_mat <- t(matrix(res[,PROT_LOG2FC], nrow=ncol(mat_TMT_MS2)))
zscores <- (l2fc_mat_TMT_MS2 - rowMeans(l2fc_mat))/rowSds(l2fc_mat)
hist(zscores, breaks=50)
table(abs(zscores) >= 3)

res[,L2FC_ZSCORE:=as.vector(zscores)]

#' 
#' # Results
#' 
DT::datatable(res[order(PROT_PADJ)])

# write it out to disk
write_tsv(res, file.path(DIR_PROC_PROT, "20180226_myResultTable_TMT_MS2.tsv"))