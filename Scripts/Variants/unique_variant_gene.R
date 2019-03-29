#'---
#' title: Create a unique variant per gene table
#' author: vyepez
#' wb:
#'  input:
#'   - var_dt: '`sm config["PROC_RESULTS"] + "/process_vcf/variant_dt.Rds"`'
#'  output:
#'   - unique_dt: '`sm config["PROC_RESULTS"] + "/process_vcf/unique_vars_gene/unique_variant_dt.Rds"`'
#'  threads: 10
#'  type: script
#'---

#+ echo=F
saveRDS(snakemake, "tmp/unique_vars.snakemake")
# snakemake <- readRDS("tmp/unique_vars.snakemake")

suppressPackageStartupMessages({
  library(data.table)
  library(BiocParallel)
})

var_dt <- readRDS(snakemake@input$var_dt)
setnames(var_dt, old = c("hgncid", "sample"), new = c("geneID", "EXOME_ID"))
var_dt <- var_dt[!is.na(geneID)]
setorder(var_dt, MAX_AF)
mt <- var_dt[, .SD[1], by = .(EXOME_ID, geneID)]
setorderv(mt, c("chr", "pos"))
saveRDS(mt, snakemake@output$unique_dt)

message('continue')

register(MulticoreParam(snakemake@threads))

bplapply(c("missense", "synonymous", "splice", "unstop", "frame-shift", "unstart", "stop", "stop_retain"), 
         function(var_type){
           print(var_type)
           sub_vt <- readRDS(paste0('/s/project/genetic_diagnosis/processed_results/process_vcf/', var_type, '_variant_dt.Rds'))
           setnames(sub_vt, old = c("hgncid", "sample"), new = c("geneID", "EXOME_ID"))
           sub_vt <- sub_vt[!is.na(geneID)]
           setorder(sub_vt, MAX_AF)
           mt <- sub_vt[, .SD[1], by = .(EXOME_ID, geneID)]
           setorderv(mt, c("chr", "pos"))
           saveRDS(mt, paste0('/s/project/genetic_diagnosis/processed_results/process_vcf/unique_vars_gene/', var_type, '_unique_variant_dt.Rds'))
})
