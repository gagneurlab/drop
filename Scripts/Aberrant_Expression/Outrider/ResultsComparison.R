#'---
#' title: Results Comparison
#' author: vyepez
#' wb:
#'  input:
#'  - results: '`sm config["PROC_RESULTS"] + "/v29_overlap/outrider/OUTRIDER_results.tsv"`'
#'  output:
#' output:
#'   html_document:
#'    code_folding: show
#'    code_download: TRUE
#'---

#+ echo=F
saveRDS(snakemake, "tmp/res_comp.snakemake")
# snakemake <- readRDS("tmp/res_comp.snakemake")

suppressPackageStartupMessages({
    library(ggplot2)
    library(ggthemes)
    library(data.table)
})

res <- fread(snakemake@input$results)
res_old <- fread("/s/project/genetic_diagnosis/processed_results/res_all_batches_th.tsv")
res_old <- res_old[LAB == 'PROKISCH']
dim(res)
dim(res_old)
res_old[HANS_CLASS == 'MITO', .N]
res[HANS_CLASS == 'MITO', .N]
# res_old[geneID %in% genes_to_check]
# res[geneID %in% genes_to_check]

new_mito_disease_genes <- setdiff(res[HANS_CLASS == 'MITO', geneID], res_old[HANS_CLASS == 'MITO', geneID])
rn <- res[geneID %in% new_mito_disease_genes, .(sampleID, geneID, l2fc, STRANDED, KNOWN_MUTATION, FIBROBLAST_ID, EXOME_ID)][order(sampleID)]
table(rn$l2fc > 0)
