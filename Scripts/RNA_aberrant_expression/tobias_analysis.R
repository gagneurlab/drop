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
setorder(res, l2fc)

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

plotExpressionRank(ods_nss, 'NDUFB9', normalized = T)
plotExpressionRank(ods_nss, 'NDUFB9', normalized = F)
plotVolcano(ods_nss, '99393')

plotExpressionRank(ods_nss, 'NDUFA10', normalized = T)
plotExpressionRank(ods_nss, 'NDUFA10', normalized = F)
plotVolcano(ods_nss, '95385')

#' ## Repeated genes
dup_genes <- res[duplicated(geneID), geneID]
dup_genes  # KIAA1586 is related to nucleic acid binding
res[geneID %in% dup_genes]  # TUEB004 and TUEB005 are siblings

plotExpressionRank(ods_nss, 'KIAA1586')
plotVolcano(ods_nss, 'TUEB004')
plotVolcano(ods_nss, 'TUEB005')

#' ## Highest ranked genes
plotExpressionRank(ods_nss, 'APOBEC3B', normalized = T)
plotExpressionRank(ods_nss, 'APOBEC3B', normalized = F)
plotVolcano(ods_nss, 'TUEB007')

plotExpressionRank(ods_nss, 'DOK5', normalized = T)
plotExpressionRank(ods_nss, 'DOK5', normalized = F)
plotVolcano(ods_nss, '97808')



#' ## Volcano plots
sat[IS_RNA_SEQ_STRANDED == TRUE, ID_Links]
plotVolcano(ods_ss, '100781R')
