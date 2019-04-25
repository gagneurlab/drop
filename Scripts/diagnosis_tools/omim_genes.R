#'---
#' title: OMIM genes
#' author: vyepez
#' wb:
#'  input:
#'   - omim_table: "/s/project/mitoMultiOmics/db_data/omim-gene-pheno-cache.RDS"
#'   - v19_dt: "/s/project/genetic_diagnosis/resource/gencode_v19_unique_gene_name.tsv"
#'   - v29_dt: "/s/project/genetic_diagnosis/resource/gencode_v29_unique_gene_name.tsv"
#'  output:
#'   - omim_genes: "/s/project/genetic_diagnosis/resource/omim_genes.Rds"
#'  type: script
#'---

#+ echo=F
saveRDS(snakemake, "tmp/omim_genes.snakemake")
# snakemake <- readRDS("tmp/omim_genes.snakemake")

suppressPackageStartupMessages({
  library(data.table)
  library(tidyr)
  library(magrittr)
  library(dplyr)
})

# Read table and subset
omim_dt = readRDS(snakemake@input$omim_table)
omim_dt = readRDS("/s/project/mitoMultiOmics/db_data/omim-gene-pheno-cache.RDS")
omim_dt <- omim_dt[PMIM != "", .(SYMBOL, PMIM, GMIM)]
omim_dt <- omim_dt[SYMBOL != ""]
omim_dt <- omim_dt[, SYMBOL := toupper(SYMBOL)]
omim_dt <- omim_dt[, .SD[1], by = SYMBOL]
dim(omim_dt)

# Add gene names from v19 and v29
v19_dt <- fread(snakemake@input$v19_dt)
v29_dt <- fread(snakemake@input$v29_dt)
v19_dt[, gene_id2 := dplyr::first(unlist(strsplit(gene_id, "\\."))), by = 1:nrow(v19_dt)]

v19_dt[, gene_name_unique := toupper(gene_name_unique)]
v29_dt[, gene_name_unique := toupper(gene_name_unique)]

# Merge annotations
omim_dt <- left_join(omim_dt, v19_dt[,.(gene_id2, gene_name_unique)], by = c("SYMBOL" = "gene_name_unique")) %>% as.data.table()
dim(omim_dt)

omim_dt <- left_join(omim_dt, v29_dt[,.(gene_id, gene_name_unique)], by = c("SYMBOL" = "gene_name_unique")) %>% as.data.table()
dim(omim_dt)

omim_dt[, c("PMIM", "GMIM") := NULL]

# Clean v29
omim_dt[, gene_v29 := SYMBOL]
omim_dt[is.na(gene_id), gene_v29 := NA]
omim_dt[is.na(gene_v29)]

omim_dt[SYMBOL == 'T', `:=` (gene_v29 = 'TBXT', gene_id = 'ENSG00000164458')]
omim_dt[SYMBOL == 'C21ORF2', `:=` (gene_v29 = 'CFAP410', gene_id = 'ENSG00000160226')]
omim_dt[SYMBOL == 'TMEM5', `:=` (gene_v29 = 'RXYLT1', gene_id = 'ENSG00000118600')]
omim_dt[SYMBOL == 'HFE2', `:=` (gene_v29 = 'HJV', gene_id = 'ENSG00000168509')]
omim_dt[SYMBOL == 'C2ORF71', `:=` (gene_v29 = 'PCARE', gene_id = 'ENSG00000179270')]
omim_dt[SYMBOL == 'C5ORF42', `:=` (gene_v29 = 'CPLANE1', gene_id = 'ENSG00000197603')]
omim_dt[SYMBOL == 'C21ORF59', `:=` (gene_v29 = 'CFAP298', gene_id = 'ENSG00000159079')]

# TCL4, PRSS2, DISC2, KCNJ18, ATXN8, MIR2861 don't exist in annotation

# Clean v19
# Merge again with new gene_id from v29
omim_dt <- left_join(omim_dt, v19_dt[,.(gene_id2, gene_name_unique)], by = c("gene_id" = "gene_id2")) %>% as.data.table
setnames(omim_dt, "gene_name_unique", "gene_v19")
omim_dt[!is.na(gene_id2), gene_v19 := SYMBOL]
omim_dt[SYMBOL == 'UQCC3', gene_v19 := 'C11ORF83']
omim_dt[SYMBOL == 'CCNQ', gene_v19 := 'FAM58A']
omim_dt[SYMBOL == 'GDF1', gene_v19 := 'CERS1']
omim_dt[SYMBOL == 'PABPN1', gene_v19 := 'BCL2L2-PABPN1']

omim_dt[, c("SYMBOL", "gene_id2") := NULL]

# Add meta info and save
omim_dt[, DISEASE := 'OMIM']
omim_dt[, N := .N, by = DISEASE]
omim_dt[, ORIGIN := 'OMIM']

omim_dt[gene_v19 != gene_v29]
omim_dt[is.na(gene_v19)]

saveRDS(omim_dt, snakemake@output$omim_genes)
