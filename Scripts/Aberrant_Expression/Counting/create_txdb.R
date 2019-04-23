#'---
#' title: Create gtf databases
#' author: mumichae, vyepez
#' wb:
#'  input:
#'   - gtex_gtf: "resources/gencode.v19.genes.patched_contigs.gtf.gz"
#'   - gencode_gtf: "/s/genomes/human/hg19/gencode29/gencode.v29lift37.annotation.gtf.gz"
#'  output:
#'  type: script
#'---

saveRDS(snakemake,  "tmp/txdb.snakemake")
# snakemake <- readRDS("tmp/txdb.snakemake")

suppressPackageStartupMessages({
    library(GenomicFeatures)
})

# v19
# downloaded from https://storage.googleapis.com/gtex_analysis_v7/reference/gencode.v19.genes.v7.patched_contigs.gtf
gtex_txdb = makeTxDbFromGFF(snakemake@input$gtex_gtf, format='gtf')
# rename chromosomes for v19
seqlevelsStyle(gtex_txdb) <- "UCSC"
#seqlevels(gtex_txdb) <- gsub("MT", "M", seqlevels(gtex_txdb))
message('v19')
message(seqlevels(gtex_txdb))
saveDb(gtex_txdb, "/s/project/genetic_diagnosis/processed_results/v19/txdb.db")

# v29
gencode_txdb = makeTxDbFromGFF(snakemake@input$gencode_gtf, format='gtf')
# Subset to include only canonical chromosomes
gencode_txdb <- keepStandardChromosomes(gencode_txdb)
message('v29')
message(seqlevels(gencode_txdb))
saveDb(gencode_txdb, "/s/project/genetic_diagnosis/processed_results/v29/txdb.db")
saveDb(gencode_txdb, "/s/project/genetic_diagnosis/processed_results/v29_overlap/txdb.db")
