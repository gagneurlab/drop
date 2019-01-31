#'---
#' title: Analysis of Tobias Haack's RNA samples
#' author: vyepez
#' wb:
#'  input: 
#'  - sample_anno: "/s/project/mitoMultiOmics/raw_data/sample_info/201812_th_sample_anno.tsv"
#'  - ods_results: "/s/project/genetic_diagnosis/processed_results/res_th.tsv"
#'  output:
#' output:
#'  html_document:
#'   code_folding: show
#'   code_download: TRUE
#'---

#+ echo=F
source('.wBuild/wBuildParser.R')
parseWBHeader("Scripts/Aberrant_Expression/tobias_analysis.R")

source("src/r/config.R")
library(OUTRIDER)
#' # Read the annotation and results tables
sat <- fread(snakemake@input[['sample_anno']])
res <- fread(snakemake@input[['ods_results']])
res <- res[LAB == "HAACK"]
setorder(res, l2fc)

ods <- readRDS("/s/project/genetic_diagnosis/processed_results/ods_kremer_th.Rds") 

#+ echo=F
plot_expected_raw_counts <- function(gene, ods){
    plot(normalizationFactors(ods[gene,]), counts(ods[gene,]), log = "xy",
         ylab = 'Raw Counts', xlab = 'Expected Counts', main = gene)
    abline(0,1)
}

#' ## How many RNAs do we have?
uniqueN(sat[ASSAY == "RNASeq", ID_Links])
#' ## How many are strand / non strand specific?
table(sat[, .(IS_RNA_SEQ_STRANDED)])
# The strand specific were counted as non-strand specific in order to merge with Kremer counts

#' ## In how many samples do we find outliers?
uniqueN(res$sampleID)

# Samples with no outliers
setdiff(sat[ASSAY == "RNASeq", ID_Links], res$sampleID)

barplot(sort(table(res$sampleID)), las = 2, ylab = "Number of expression outliers", cex.lab = 1.2); grid()

#' ## Results table
DT::datatable(res, caption = "Expression outlier results table", style = 'bootstrap')

#' ## Mito disease or Mitocarta genes
DT::datatable(res[!is.na(MITOGENE_CATEGORY) | MITOCARTA == TRUE])
# Three promising genes: NDUFA10, NDUFA2 and NDUFB9

plotExpressionRank(ods, 'NDUFB9', normalized = T)
plot_expected_raw_counts("NDUFB9", ods)
plotVolcano(ods, '99393')

plotExpressionRank(ods, 'NDUFA10', normalized = T)
plot_expected_raw_counts("NDUFA10", ods)
plotVolcano(ods, '95385')

plotExpressionRank(ods, 'NDUFA2', normalized = T)

#' ## Repeated genes
dup_genes <- res[duplicated(geneID), geneID]
dup_genes  # KIAA1586 is related to nucleic acid binding
res[geneID %in% dup_genes]  # TUEB004 and TUEB005 are siblings

plotExpressionRank(ods, 'KIAA1586')
plotVolcano(ods, 'TUEB004')
plotVolcano(ods, 'TUEB005')



