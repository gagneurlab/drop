#'---
#' title: Check possibly contaminated RNAs
#' author: vyepez
#' wb:
#'  input:
#'   - kremer_ods: '/s/project/genetic_diagnosis/processed_results/ods_kremer_counts.Rds'
#'  output:
#' output: 
#'   html_document:
#'    code_folding: show
#'    code_download: TRUE
#'---

#' ### Robert sent a list of the first RNAs sequenced that were done under different protocols and might have a negative impact in the analysis. The goal is to decide if they should be removed or not.

library(OUTRIDER)
library(magrittr)
library(dplyr)
sa <- fread("../sample_annotation/Data/sample_annotation.tsv")
bad_rnas <- c(35834,57415,61695,61982,65937,66623,69245,69248,69456,69878,70038,70041,70476,70477,70478,70479,72748,74123,74172)
# All are fibros except for 70476 and 70477. We do not have RNA 70477 and 69878.
length(bad_rnas)
sa[RNA_ID %in% bad_rnas, .(FIBROBLAST_ID, EXOME_ID, RNA_ID, GENOME_ID, KNOWN_MUTATION, TISSUE, DISEASE, MAX_RNA_TO_EXOME_IDENTITY, COMMENT, RNA_PERSON)][order(RNA_ID)]


ods <- readRDS(snakemake@input[['kremer_ods']])
bad_rnas %in% colnames(ods) %>% sum
bad_rnas[! bad_rnas %in% colnames(ods)]
sa[RNA_ID %in% bad_rnas[! bad_rnas %in% colnames(ods)]]  # Corresponds to the blood and transduced

colData(ods)$bad_rna = colData(ods)$sampleID %in% bad_rnas

#' ## Check if they cluster before and after normalizing
#+ fig.width=12, fig.height=12
plotCountCorHeatmap(ods, normalized = F, nCluster = 2, rowCoFactor = "bad_rna")
# ods <- filterExpression(ods, '/s/genomes/human/hg19/ucsc/ucsc.translated.gtf')

#+ fig.width=12, fig.height=12
plotCountCorHeatmap(ods, normalized = T, nCluster = 2, rowCoFactor = "bad_rna")
# The samples that cluster are NHDFs

#' ## Check for excessive number of outlier genes
res <- OUTRIDER::results(ods)
tab = table(res$sampleID) %>% sort

#+ fig.width=14
barplot(tab, col = names(tab)%in% bad_rnas, las = 2, ylab = 'Outlier genes')

#' ### Conclusion: don't remove the samples.
