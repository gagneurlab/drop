#'---
#' title: Create GeneID-GeneName mapping
#' author: mumichae
#' wb:
#'  params:
#'   - tmpdir: '`sm drop.getMethodPath(METHOD, "tmp_dir")`'
#'  input:
#'   - gtf: '`sm lambda wildcards: parser.getGeneAnnotationFile(wildcards.annotation) `'
#'  output:
#'   - gene_name_mapping: '`sm parser.getProcDataDir() + "/mae/gene_name_mapping_{annotation}.tsv"`'
#'  type: script
#'---

saveRDS(snakemake,  paste0(snakemake@config$tmpdir, "/MAE/gene_map.snakemake"))
# snakemake <- readRDS(paste0(snakemake@config$tmpdir, "/MAE/gene_map.snakemake"))

suppressPackageStartupMessages({
  library(rtracklayer)
  library(data.table)
  library(magrittr)
  library(tidyr)
})

gtf_dt <- import(snakemake@input$gtf) %>% as.data.table
if (!"gene_name" %in% colnames(gtf_dt)) {
  gtf_dt[gene_name := gene_id]
}
gtf_dt <- gtf_dt[type == "gene", .(seqnames, start, end, strand, gene_id, gene_name, gene_type)]

# make gene_names unique
gtf_dt[, N := 1:.N, by = gene_name] # warning message
gtf_dt[, gene_name_orig := gene_name]
gtf_dt[N > 1, gene_name := paste(gene_name, N, sep = '_')]
gtf_dt[, N := NULL]

fwrite(gtf_dt, snakemake@output$gene_name_mapping, na = NA)
