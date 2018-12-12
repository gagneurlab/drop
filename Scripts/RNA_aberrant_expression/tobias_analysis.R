#'---
#' title: Analysis of Tobias Haack's RNA samples
#' author: vyepez
#' wb:
#'  input: 
#'  - sample_anno: "/s/project/mitoMultiOmics/raw_data/sample_info/201812_th_sample_anno.tsv"
#'  - ods_results: "/s/project/genetic_diagnosis/processed_results/res_all_batches_th.tsv"
#'  output:
#' output:
#'  html_document:
#'   code_folding: show
#'   code_download: TRUE
#'---

#+ echo=F
source('.wBuild/wBuildParser.R')
parseWBHeader("Scripts/RNA_aberrant_expression/tobias_analysis.R")

source("src/r/config.R")
library(OUTRIDER)

#' # Read the annotation and results tables
sat <- fread(snakemake@input[['sample_anno']])

res <- fread(snakemake@input[['ods_results']])
res <- res[LAB == "HAACK"]

ods_ss <- readRDS("/s/project/genetic_diagnosis/processed_results/ods_batches2_3_4_th_ss.Rds")
ods_nss <- readRDS("/s/project/genetic_diagnosis/processed_results/ods_batches0_1_th_nss.Rds")

#' ## How many RNAs do we have?
uniqueN(sat[ASSAY == "RNASeq", ID_Links])
#' ## How many are strand / non strand specific?
table(sat[, .(IS_RNA_SEQ_STRANDED)])


#' ## How many outliers do we find?
uniqueN(res$sampleID)

# Samples with no outliers
setdiff(sat[ASSAY == "RNASeq", ID_Links], res$sampleID)

barplot(sort(table(res$sampleID)), las = 2, ylab = "Number of expression outliers", cex.lab = 1.2); grid()

#' ## Results table
DT::datatable(res, caption = "Expression outlier results table", style = 'bootstrap')

#' ## Mito disease or Mitocarta genes
DT::datatable(res[!is.na(MITOGENE_CATEGORY) | MITOCARTA == TRUE])

plotExpressionRank(ods_nss, 'ENSG00000203667.5', main = "COX20")

plotExpressionRank(ods_nss, 'ENSG00000130414.7', normalized = F)

plotExpressionRank(ods_nss, 'ENSG00000130414.7', normalized = F)

#' ## Repeated genes
dup_genes <- res[duplicated(gene_name), gene_name]
dup_genes  # KIAA1586 is related to nucleic acid binding
res[gene_name %in% dup_genes]  # TUEB004 and TUEB005 are siblings

plotExpressionRank(ods_nss, 'ENSG00000168116.9', main = "KIAA1586")


#' ## Volcano plots
sat[IS_RNA_SEQ_STRANDED == TRUE, ID_Links]
plotVolcano(ods_ss, '100781R')
