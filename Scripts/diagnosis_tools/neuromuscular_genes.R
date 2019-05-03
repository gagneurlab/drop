#'---
#' title: Neuromuscular disease genes
#' author: vyepez
#' wb:
#'  input:
#'   - v19_dt: "/s/project/genetic_diagnosis/resource/gencode_v19_unique_gene_name.tsv"
#'   - v29_dt: "/s/project/genetic_diagnosis/resource/gencode_v29_unique_gene_name.tsv"
#'  output:
#'   - neuromuscular_genes: "/s/project/genetic_diagnosis/resource/neuromuscular_genes.tsv"
#'  type: script
#'---

#+ echo=F
saveRDS(snakemake, "tmp/neuromuscular_genes.snakemake")
# snakemake <- readRDS("tmp/neuromuscular_genes.snakemake")

suppressPackageStartupMessages({
  library(data.table)
  library(tidyr)
  library(magrittr)
  library(dplyr)
})

DIR_ext_genes <- "/s/project/mitoMultiOmics/raw_data/gene_info/external"

# Neuromuscular Genes
neuro_genes_dt <- read.table(file.path(DIR_ext_genes, "neuromuscular_genes_gonorazky.txt"),
                             skip = 1, col.names = 'gene_id') %>% as.data.table()

neuro_genes_dt[, DISEASE := 'Neuromuscular']

# Add gene names from v19 and v29
v19_dt <- fread(snakemake@input$v19_dt)
v29_dt <- fread(snakemake@input$v29_dt)
# v19_dt <- fread("/s/project/genetic_diagnosis/resource/gencode_v19_unique_gene_name.tsv")
# v29_dt <- fread("/s/project/genetic_diagnosis/resource/gencode_v29_unique_gene_name.tsv")
v19_dt[, gene_id2 := dplyr::first(unlist(strsplit(gene_id, "\\."))), by = 1:nrow(v19_dt)]

dim(neuro_genes_dt)
neuro_genes_dt <- left_join(neuro_genes_dt, v19_dt[, .(gene_id2, gene_name)], 
                              by = c("gene_id" = "gene_id2")) %>% as.data.table()
dim(neuro_genes_dt)

setnames(neuro_genes_dt, "gene_name", "gene_v19")
neuro_genes_dt[gene_id == 'ENSG00000100836', gene_v19 := 'BCL2L2-PABPN1']

neuro_genes_dt <- left_join(neuro_genes_dt, v29_dt[, .(gene_id, gene_name)], by = "gene_id") %>% as.data.table()
dim(neuro_genes_dt)
setnames(neuro_genes_dt, "gene_name", "gene_v29")
neuro_genes_dt[gene_id == 'ENSG00000182500', `:=` (gene_id = 'ENSG00000276045',
                                                   gene_v29 = 'ORAI1')]  # gene id changed

neuro_genes_dt[, N := .N, by = DISEASE]
neuro_genes_dt[, ORIGIN := 'Gonorazky']

neuro_genes_dt[gene_v19 != gene_v29]

fwrite(neuro_genes_dt, snakemake@output$neuromuscular_genes)

