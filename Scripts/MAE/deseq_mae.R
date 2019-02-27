#'---
#' title: Get MAE results
#' author: vyepez
#' wb:
#'  input:
#'   - mae_counts: '`sm config["PROC_DATA"] + "/mae/{vcf}-{rna}.Rds"`'
#'  output:
#'   - mae_res: '`sm config["PROC_RESULTS"] + "/mae/samples/{vcf}-{rna}_res.Rds"`'
#'  type: script
#'---

saveRDS(snakemake, 'tmp/res_mae.Rds')
# snakemake <- readRDS(snakemake, 'tmp/res_mae.Rds')

suppressPackageStartupMessages({
    devtools::load_all("../mae/")
})

mae_raw <- readRDS(snakemake@input$mae_counts)

# Function from MAE own pkg
rmae <- run_deseq_all_mae(mae_raw)

saveRDS(rmae, snakemake@output$mae_res)
