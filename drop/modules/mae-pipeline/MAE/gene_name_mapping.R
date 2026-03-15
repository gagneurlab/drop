#'---
#' title: Create GeneID-GeneName mapping
#' author: mumichae
#' wb:
#'  log:
#'   - snakemake: '`sm str(tmp_dir / "MAE" / "{annotation}.Rds")`'
#'  input:
#'   - gtf: '`sm lambda w: cfg.genome.getGeneAnnotationFile(w.annotation) `'
#'  output:
#'   - gene_name_mapping: '`sm cfg.getProcessedDataDir() + "/mae/gene_name_mapping_{annotation}.tsv"`'
#'  type: script
#'---

saveRDS(snakemake, snakemake@log$snakemake)

suppressPackageStartupMessages({
  library(rtracklayer)
  library(data.table)
  library(magrittr)
  library(tidyr)
})

gtf_dt <- import(snakemake@input$gtf) %>% as.data.table
if (!"gene_name" %in% colnames(gtf_dt)) {
  gtf_dt[, gene_name := gene_id]
}
if('gene_biotype' %in% colnames(gtf_dt))
   setnames(gtf_dt, 'gene_biotype', 'gene_type')
gtf_dt <- gtf_dt[type == "gene", .(seqnames, start, end, strand, gene_id, gene_name, gene_type)]

# make gene_names unique, keeping the originals for reference.
# make.unique avoids collisions that occur with the N-counter approach when
# a gene is already naturally named e.g. "GENEX_2" (two different genes would
# otherwise both map to the same name).
gtf_dt[, gene_name_orig := gene_name]
gtf_dt[, gene_name := make.unique(gene_name, sep = '_')]

fwrite(gtf_dt, snakemake@output$gene_name_mapping, na = NA)
