#'---
#' title: OUTRIDER Results
#' author: mumichae
#' wb:
#'  input:
#'   - ods: '`sm config["PROC_RESULTS"] + "/{annotation}/outrider/{dataset}/ods.Rds"`'
#'   - functions: "Scripts/_functions/gene_annotation/add_gene_info_cols.R"
#'  output:
#'   - results: '`sm config["PROC_RESULTS"] + "/{annotation}/outrider/{dataset}/OUTRIDER_results.tsv"`'
#'   - results_all: '`sm config["PROC_RESULTS"] + "/{annotation}/outrider/{dataset}/OUTRIDER_results_all.Rds"`'
#'   - results_public: "/s/public_webshare/project/genetic_diagnosis/results/{annotation}/OUTRIDER_results_{dataset}.tsv"
#'  type: script
#'---

#+ echo=F
saveRDS(snakemake, "tmp/outrider_results.snakemake")
# snakemake <- readRDS("tmp/outrider_results.snakemake")

suppressPackageStartupMessages({
    library(OUTRIDER)
    library(SummarizedExperiment)
    library(ggplot2)
    library(data.table)
    library(dplyr)
})

source(snakemake@input$functions)

ods <- readRDS(snakemake@input$ods)
res <- OUTRIDER::results(ods, all = TRUE)
res[, FC := round(2^l2fc, 2)]
res[, geneID := toupper(geneID)]
res <- add_gene_type(res, gene_name_col = 'geneID')
saveRDS(res[,.(geneID, sampleID, pValue, padjust, zScore, l2fc, rawcounts, normcounts, meanCorrected, theta, aberrant, AberrantBySample, AberrantByGene, padj_rank, FC, gene_type)], snakemake@output$results_all)

# Subset to significant results
res <- res[padjust <= .05]
res <- add_all_gene_info(res, gene_name_col = 'geneID', dis_genes = F, gene_type = F)  # gene_type already added before

# Add sample annotation
sa <- fread(snakemake@config$SAMPLE_ANNOTATION)
res <- left_join(res, sa[, .(RNA_ID, FIBROBLAST_ID, EXOME_ID, PEDIGREE, KNOWN_MUTATION,
                             CANDIDATE_GENE, BATCH, COMMENT, PROTEOME_ID, DISEASE, RNA_PERSON)],
                 by = c("sampleID" = "RNA_ID")) %>% as.data.table

res[, tp_sample := as.character(any(geneID == KNOWN_MUTATION)), by = sampleID]
res[is.na(KNOWN_MUTATION), tp_sample := "Unsolved"]

# Save results 
fwrite(res, snakemake@output$results, sep = "\t", quote = F)
fwrite(res, snakemake@output$results_public, sep = "\t", quote = F)
