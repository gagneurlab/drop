#'---
#' title: Mito disease genes
#' author: vyepez
#' wb:
#'  input:
#'   - hans_table: "/s/project/mitoMultiOmics/raw_data/gene_info/mitochondrial_disorder_genes_prokisch_mayr_cleaned.tsv"
#'   - v19_dt: "/s/project/genetic_diagnosis/resource/gencode_v19_unique_gene_name.tsv"
#'   - v29_dt: "/s/project/genetic_diagnosis/resource/gencode_v29_unique_gene_name.tsv"
#'  output:
#'   - mito_genes: "/s/project/genetic_diagnosis/resource/mito_disease_genes.tsv"
#'  type: script
#'---

#+ echo=F
saveRDS(snakemake, "tmp/md_genes.snakemake")
# snakemake <- readRDS("tmp/md_genes.snakemake")

suppressPackageStartupMessages({
  library(data.table)
  library(tidyr)
  library(magrittr)
  library(dplyr)
})

# Read table
hans_dt <- fread(snakemake@input$hans_table)
mito_genes_dt <- hans_dt[DISEASE == 'MITO', .(HGNC_GENE_NAME)]
dim(mito_genes_dt)

# Add gene names from v19 and v29
v19_dt <- fread(snakemake@input$v19_dt)
v29_dt <- fread(snakemake@input$v29_dt)
# v19_dt <- fread("/s/project/genetic_diagnosis/resource/gencode_v19_unique_gene_name.tsv")
# v29_dt <- fread("/s/project/genetic_diagnosis/resource/gencode_v29_unique_gene_name.tsv")
v19_dt[, gene_id2 := dplyr::first(unlist(strsplit(gene_id, "\\."))), by = 1:nrow(v19_dt)]

v19_dt[, gene_name := toupper(gene_name)]
v29_dt[, gene_name_unique := toupper(gene_name_unique)]

mito_genes_dt <- left_join(mito_genes_dt, v19_dt[,.(gene_id2, gene_name)], by = c("HGNC_GENE_NAME" = "gene_name")) %>% as.data.table()
dim(mito_genes_dt)

mito_genes_dt <- left_join(mito_genes_dt, v29_dt[,.(gene_id, gene_name_unique)], by = c("HGNC_GENE_NAME" = "gene_name_unique")) %>% as.data.table()
dim(mito_genes_dt)

mito_genes_dt[is.na(gene_id)]
mito_genes_dt[, gene_v29 := HGNC_GENE_NAME]
mito_genes_dt[HGNC_GENE_NAME == 'FDX1L', `:=` (gene_v29 = 'FDX2', gene_id = 'ENSG00000267673')]
mito_genes_dt[HGNC_GENE_NAME == 'ATP5MD_USMG5', `:=` (gene_v29 = 'ATP5MD', gene_id = 'ENSG00000173915')]

mito_genes_dt[is.na(gene_id2)]
mito_genes_dt[, gene_v19 := HGNC_GENE_NAME]
mito_genes_dt[is.na(gene_id2), gene_v19 := NA]

mito_genes_dt[HGNC_GENE_NAME == 'COQ8A', gene_v19 := 'ADCK3']
mito_genes_dt[HGNC_GENE_NAME == 'COQ8B', gene_v19 := 'ADCK4']
mito_genes_dt[HGNC_GENE_NAME == 'ATP5F1E', gene_v19 := 'ATP5E']
mito_genes_dt[HGNC_GENE_NAME == 'ATP5F1A', gene_v19 := 'ATP5A1']
mito_genes_dt[HGNC_GENE_NAME == 'ATP5F1D', gene_v19 := 'ATP5D']
mito_genes_dt[HGNC_GENE_NAME == 'FDX1L', gene_v19 := 'FDX2']
mito_genes_dt[HGNC_GENE_NAME == 'TWNK', gene_v19 := 'C10ORF2']
mito_genes_dt[HGNC_GENE_NAME == 'ATP5MD_USMG5', gene_v19 := 'ATP5MD']
mito_genes_dt[HGNC_GENE_NAME == 'UQCC3', gene_v19 := 'C11ORF83']
mito_genes_dt[HGNC_GENE_NAME == 'COA7', gene_v19 := 'SELRC1']
mito_genes_dt[HGNC_GENE_NAME == 'NAXE', gene_v19 := 'APOA1BP']
mito_genes_dt[HGNC_GENE_NAME == 'NAXD', gene_v19 := 'CARKD']
mito_genes_dt[HGNC_GENE_NAME == 'GATB', gene_v19 := 'PET112']
mito_genes_dt[HGNC_GENE_NAME == 'NDUFAF8', gene_v19 := 'C17ORF89']
mito_genes_dt[HGNC_GENE_NAME == 'MRM2', gene_v19 := 'FTSJ2']

mito_genes_dt[, gene_id2 := NULL]
mito_genes_dt[, HGNC_GENE_NAME := NULL]

mito_genes_dt[, DISEASE := 'MITO']
mito_genes_dt[, N := .N, by = DISEASE]
mito_genes_dt[, ORIGIN := 'Hans Mayr']

mito_genes_dt[gene_v19 != gene_v29]
mito_genes_dt[is.na(gene_v19)]

fwrite(mito_genes_dt, snakemake@output$mito_genes)
