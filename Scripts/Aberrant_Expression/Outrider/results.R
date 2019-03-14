#'---
#' title: OUTRIDER Results
#' author: mumichae
#' wb:
#'  input:
#'   - ods: '`sm config["PROC_RESULTS"] + "/{annotation}/outrider/{dataset}/ods.Rds"`'
#'   - functions: "Scripts/_functions/gene_annotation/add_gene_info_cols.R"
#'  output:
#'   - results: '`sm config["PROC_RESULTS"] + "/{annotation}/outrider/{dataset}/OUTRIDER_results.tsv"`'
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
res <- results(ods)
res[, FC := 2^l2fc]
res <- add_all_gene_info(res, gene_name = 'geneID', dis_genes = F)

# Add sample annotation
sa <- fread(snakemake@config$SAMPLE_ANNOTATION)
res[, geneID := toupper(geneID)]
res <- left_join(res, sa[, .(RNA_ID, FIBROBLAST_ID, EXOME_ID, PEDIGREE, KNOWN_MUTATION,
                             CANDIDATE_GENE, BATCH)],
                 by = c("sampleID" = "RNA_ID")) %>% as.data.table

res[, TP := as.character(geneID == KNOWN_MUTATION)]
res[is.na(KNOWN_MUTATION), TP := "Unsolved"]

res[, tp_sample := as.character(any(geneID == KNOWN_MUTATION)), by = sampleID]
res[is.na(KNOWN_MUTATION), tp_sample := "Unsolved"]

fwrite(res, snakemake@output$results, sep = "\t", quote = F)
fwrite(res, snakemake@output$results_public, sep = "\t", quote = F)
