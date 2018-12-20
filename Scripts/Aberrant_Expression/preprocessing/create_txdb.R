#'---
#' title: Create gtf databases
#' author: mumichae, vyepez
#' wb:
#'  input:
#'   - gtex_gtf: "resources/gencode.v19.genes.patched_contigs.gtf.gz"
#'   - gencode_gtf: "/s/genomes/human/hg19/gencode29/gencode.v29lift37.annotation.gtf.gz"
#'  output:
#'   - gtex_txdb: "resources/gencode.v19.genes.patched_contigs.Db"
#'   - gencode_txdb: "/s/genomes/human/hg19/gencode29/gencode.v29lift37.annotation.Db"
#'  type: script
#'---

saveRDS(snakemake,  "tmp/txdb.snakemake")
suppressPackageStartupMessages({
    library(GenomicFeatures)
})

# v19
# downloaded from https://storage.googleapis.com/gtex_analysis_v7/reference/gencode.v19.genes.v7.patched_contigs.gtf
gtex_txdb = makeTxDbFromGFF(snakemake@input$gtex_gtf, format='gtf')
# rename chromosomes for v19
seqlevels(gtex_txdb) <- paste0("chr", seqlevels(gtex_txdb))
seqlevels(gtex_txdb) <- gsub("MT", "M", seqlevels(gtex_txdb))
saveDb(gtex_txdb, snakemake@output$gtex_txdb)

# v29
gencode_txdb = makeTxDbFromGFF(snakemake@input$gencode_gtf, format='gtf')
# Subset to include only canonical chromosomes
std_chr = paste0('chr',c('X','Y','M',1:22))
seqlevels(gencode_txdb) <- std_chr  # Subset the whole txdb object
saveDb(gencode_txdb, snakemake@output$gencode_txdb)
