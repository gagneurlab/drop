#'---
#' title: Create count annotation objects
#' author: mumichae, vyepez
#' wb:
#'  input:
#'   - gtex_gene_mapping: "resources/GENCODEv19Mapping.RDS"
#'   - gtex_txdb: "resources/gencode.v19.genes.patched_contigs.Db"
#'   - gencode_gtf: "/s/genomes/human/hg19/gencode29/gencode.v29lift37.annotation.gtf.gz"
#'  output:
#'   - gtex_dt: "resources/gencode_v19_unique_gene_name.tsv"
#'   - gencode_dt: "resources/gencode_v29_unique_gene_name.tsv"
#'  type: script
#'---

saveRDS(snakemake,  "tmp/unique_name.snakemake")
suppressPackageStartupMessages({
    library(BSgenome.Hsapiens.UCSC.hg19)
    library(GenomicFeatures)
    library(GenomicRanges)
    library(data.table)
    library(magrittr)
    library(dplyr)
    library(tidyr)
})


# v19
gtex_txdb <- loadDb(snakemake@input$gtex_txdb)
genes_en = sort(genes(gtex_txdb))
gene_mapping <- readRDS(snakemake@input$gtex_gene_mapping)

genes_dt = as.data.table(genes_en)
genes_dt = merge(genes_dt, gene_mapping, by = "gene_id")

## We might need to add _2, _3, ... for repeated gene names
# We can use the following code
genes_dt[, N := 1:.N, by = gene_name]
genes_dt[, gene_new_name := gene_name]
genes_dt[N > 1, gene_new_name := paste(gene_name, N, sep = "_")]
genes_dt[, gene_name := gene_new_name]
genes_dt[, N := NULL]
genes_dt[, gene_new_name := NULL]
genes_dt[, gene_id_unique := gene_id]
fwrite(genes_dt, snakemake@output$gtex_dt)


# v29
gtf_or <- rtracklayer::import(snakemake@input$gencode_gtf) %>% as.data.table
gtf_dt <- copy(gtf_or)
gtf_dt <- gtf_dt[type == "gene", .(seqnames, start, end, strand, gene_id, gene_name, gene_type, gene_status)]
gtf_dt <- gtf_dt[seqnames %in% GenomeInfoDb::standardChromosomes(BSgenome.Hsapiens.UCSC.hg19)]
setnames(gtf_dt, "gene_id", "gene_id_unique")
gtf_dt <- separate(gtf_dt, "gene_id_unique", into = "gene_id", sep = "\\.", remove = F)

dup_genes <- gtf_dt[duplicated(gtf_dt$gene_name), gene_name] # Get genes that appear at least twice
# Get genes that appear more than twice
repeated_genes <- names(table(gtf_dt[gene_name %in% dup_genes, gene_name])[table(gtf_dt[gene_name %in% dup_genes, gene_name]) > 1])

# rename duplicate gene names
gtf_dt[, N := 1:.N, by = gene_name] # warning message
gtf_dt[, gene_name_unique := gene_name]
gtf_dt[N > 1, gene_name_unique := paste(gene_name, N, sep = '_')]
gtf_dt[, N := NULL]

# check if successful
gtf_dt[gene_name %in% repeated_genes]
gtf_dt[duplicated(gtf_dt$gene_id)]  # 41 X-Y paralog genes
dup_genes <- gtf_dt[duplicated(gtf_dt$gene_name), gene_name]
# View(gtf_dt[gene_name %in% dup_genes])
fwrite(gtf_dt, snakemake@output$gencode_dt, sep = '\t')

