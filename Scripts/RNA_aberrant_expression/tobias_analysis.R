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

# 1. commit tobias_analysis.R from laptop
# 2. pull from serve
# 3. in outrider.R change the gene_id to gene_name
# 4. save ods_ss on top
# 5. recompute results and save them
# 6. check plots work
# 7. run snakemake

source("src/r/config.R")
library(OUTRIDER)
DIR_lrz = "../../../LRZ Sync+Share/LMU/TH"
#' # Read the annotation and results tables
sat <- fread(snakemake@input[['sample_anno']])
sat <- fread(file.path(DIR_lrz, "201812_th_sample_anno.tsv"))
res <- fread(snakemake@input[['ods_results']])
res <- fread(file.path(DIR_lrz, "res_all_batches_th.tsv"))
res <- res[LAB == "HAACK"]
setorder(res, l2fc)

ods_ss <- readRDS("/s/project/genetic_diagnosis/processed_results/ods_batches2_3_4_th_ss.Rds")
ods_ss <- readRDS(file.path(DIR_lrz, "ods_batches2_3_4_th_ss.Rds"))
ods_nss <- readRDS("/s/project/genetic_diagnosis/processed_results/ods_batches0_1_th_nss.Rds")
ods_nss <- readRDS(file.path(DIR_lrz, "ods_batches0_1_th_nss.Rds"))

#'+echo=F
plot_expected_raw_counts <- function(gene, ods){
    plot(normalizationFactors(ods[gene,]), counts(ods[gene,]), log = "xy",
         ylab = 'Raw Counts', xlab = 'Expected Counts', main = gene)
    abline(0,1)
}

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
plot_expected_raw_counts("NDUFB9", ods_nss)
plotVolcano(ods_nss, '99393')

plotExpressionRank(ods_nss, 'NDUFA10', normalized = T)
plotExpressionRank(ods_nss, 'NDUFA10', normalized = F)

plot_expected_raw_counts("NDUFA10", ods_nss)
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
plot_expected_raw_counts('APOBEC3B', ods_nss)
plotVolcano(ods_nss, 'TUEB007')

plotExpressionRank(ods_nss, 'DOK5', normalized = T)
plotExpressionRank(ods_nss, 'DOK5', normalized = F)
plot_expected_raw_counts('DOK5', ods_nss)
plotVolcano(ods_nss, '97808')

plotExpressionRank(ods_ss, 'ENSG00000119321.4', normalized = T, main = 'FKBP15')
plotExpressionRank(ods_ss, 'ENSG00000119321.4', normalized = F, main = 'FKBP15')
plot_expected_raw_counts('ENSG00000119321.4', ods_ss)

plotExpressionRank(ods_nss, 'AAGAB', normalized = T)
plotExpressionRank(ods_nss, 'AAGAB', normalized = F)
plot_expected_raw_counts('AAGAB', ods_nss)
plotVolcano(ods_nss, '93067')

plotExpressionRank(ods_nss, 'TSSC4', normalized = T)
plotExpressionRank(ods_nss, 'TSSC4', normalized = F)
plot_expected_raw_counts('TSSC4', ods_nss)
plotVolcano(ods_nss, '96169')

plotExpressionRank(ods_nss, 'ACSF3', normalized = T)
plotExpressionRank(ods_nss, 'ACSF3', normalized = F)
plot_expected_raw_counts('ACSF3', ods_nss)
plotVolcano(ods_nss, '97758')
