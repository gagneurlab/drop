#'---
#' title: MAE Results table
#' author: vyepez
#' wb:
#'  input:
#'   - mae_res: '`sm expand(config["PROC_RESULTS"] + "/mae/samples/{id}_res.Rds", id = config["mae_ids"])`'
#'  output:
#'   - res_signif_all: '`sm config["PROC_RESULTS"] + "/mae/MAE_results.Rds"`'
#'   - res_signif_rare: '`sm config["PROC_RESULTS"] + "/mae/MAE_results_rare.Rds"`'
#' output: 
#'   html_document:
#'    code_folding: show
#'    code_download: TRUE
#'---

#+ echo=F
saveRDS(snakemake, 'tmp/mae_res_all.Rds')
# snakemake <- readRDS('tmp/mae_res_all.Rds')
suppressPackageStartupMessages({
    library(data.table)
    library(magrittr)
    library(ggplot2)
    library(tidyr)
})

source("Scripts/_functions/gene_annotation/add_gene_info_cols.R")

#' ## Read all mae files
# maes <- list.files("/s/project/genetic_diagnosis/processed_results/mae/samples/", full.names = TRUE)

res <- lapply(snakemake@input$mae_res, function(m){
    rt <- readRDS(m)
    rt <- rt[padj < .05 & alt_freq > .8]
    return(rt)
}) %>% rbindlist()


res <- separate(res, 'sample', into = c('EXOME_ID', 'RNA_ID'), sep = "-", remove = FALSE)

# Add gene info
setnames(res, "hgncid", "gene_name")
res <- add_all_gene_info(res, dis_genes = F)

# Add sample annotation
sa <- fread(snakemake@config$SAMPLE_ANNOTATION)
res[, gene_name := toupper(gene_name)]
res <- left_join(res, sa[, .(RNA_ID, FIBROBLAST_ID, EXOME_ID, PEDIGREE, KNOWN_MUTATION,
                             CANDIDATE_GENE, BATCH, COMMENT)],
                 by = c("EXOME_ID", "RNA_ID")) %>% as.data.table

saveRDS(res, snakemake@output$res_signif_all)

res_rare <- res[MAX_AF < .001 | is.na(MAX_AF)]
saveRDS(res_rare, snakemake@output$res_signif_rare)

#' ### Download results tables
write.table(res_rare, "/s/public_webshare/project/genetic_diagnosis/results/MAE_results_rare.tsv", sep = "\t", quote = F, row.names = F)

#' [Download MAE rare results table](https://i12g-gagneurweb.informatik.tu-muenchen.de/project/genetic_diagnosis/results/MAE_results_rare.tsv)
DT::datatable(res_rare, caption = "MAE results", style = 'bootstrap')

#' ## Plot
hist(res$alt_freq, breaks = 20)
hist(res_rare$alt_freq, breaks = 20)

