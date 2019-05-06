#'---
#' title: Create count Granges from annotation
#' author: mumichae, vyepez
#' wb:
#'  input:
#'   - gtex_txdb: '`sm config["PROC_RESULTS"] + "/v19/txdb.db"`'
#'   - gencode_txdb: '`sm config["PROC_RESULTS"] + "/v29/txdb.db"`'
#'  output:
#'   - gtex_op: '`sm config["PROC_RESULTS"] + "/v19/counts/exons_by_gene_op.Rds"`'
#'   - gencode_op: '`sm config["PROC_RESULTS"] + "/v29/counts/exons_by_gene_op.Rds"`'
#'   - gencode_ov_op: '`sm config["PROC_RESULTS"] + "/v29_overlap/counts/exons_by_gene_op.Rds"`'
#'  type: script
#'---

#TODO: change name of output files for later input! rds vs. Rds

saveRDS(snakemake,  "tmp/count_ranges.snakemake")
# snakemake <- readRDS("tmp/count_ranges.snakemake")
suppressPackageStartupMessages({
    library(GenomicFeatures)
    library(GenomicRanges)
})

gtex_txdb <- loadDb(snakemake@input$gtex_txdb)
seqlevelsStyle(gtex_txdb) <- "UCSC"
gencode_txdb <- loadDb(snakemake@input$gencode_txdb)
gencode_txdb <- keepStandardChromosomes(gencode_txdb)

invert_strand <- function(exons) {
    exons_op <- data.table::copy(exons)
    if (class(exons) == "GRanges") {
        strand(exons_op) <- ifelse(strand(exons) == '+', '-', '+')
    } else if (class(exons) %in% c("GRangesList", "CompressedGRangesList")) {
        unlisted <- unlist(exons_op) # create genomic ranges object
        strand(unlisted) <- ifelse(strand(unlisted) == '+', '-', '+')
        exons_op <- relist(unlisted, exons_op) # change back to GRList
    } else {
        message(class(exons))
    }
    exons_op
}

# gtex gencode 19
gtex_op <- invert_strand(exonsBy(gtex_txdb, by = "gene"))
saveRDS(gtex_op, snakemake@output$gtex_op)

# gencode 29
gencode_op <- invert_strand(exonsBy(gencode_txdb, by = "gene"))
saveRDS(gencode_op, snakemake@output$gencode_op)
saveRDS(gencode_op, snakemake@output$gencode_ov_op)

