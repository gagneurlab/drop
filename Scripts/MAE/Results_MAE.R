#'---
#' title: MAE Results table
#' author: vyepez
#' wb:
#'  input:
#'   - mae_res: '`sm expand(config["PROC_RESULTS"] + "/mae/samples/{id}_res.Rds", id = config["mae_ids"])`'
#'  output:
#'   - res_signif_all: '`sm config["PROC_RESULTS"] + "/mae/MAE_results.Rds"`'
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
})

#' ## Read all mae files
# maes <- list.files("/s/project/genetic_diagnosis/processed_results/mae/samples/", full.names = TRUE)

res_all <- lapply(snakemake@input$mae_res, function(m){
    rt <- readRDS(m)
    rt <- rt[padj < .05 & alt_freq >= .8]
    return(rt)
}) %>% rbindlist()

saveRDS(res_all, snakemake@output$res_signif_all)


#' ### Download results table
write.table(res_all, "/s/public_webshare/project/genetic_diagnosis/results/MAE_results.tsv", sep = "\t", quote = F, row.names = F)

#' [Download MAE results table](https://i12g-gagneurweb.informatik.tu-muenchen.de/project/genetic_diagnosis/results/MAE_results.tsv)
# DT::datatable(res_all, caption = "MAE results", style = 'bootstrap')

#' ## Plot
hist(res_all$alt_freq, breaks = 20)

