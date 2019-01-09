#'---
#' title: Create count Granges from annotation
#' author: mumichae, vyepez
#' wb:
#'  input:
#'   - gtex_txdb: '`sm config["PROC_RESULTS"] + "/counts/v19/txdb.Rds"`'
#'   - gencode_txdb: '`sm config["PROC_RESULTS"] + "/counts/v29/txdb.Rds"`'
#'  output:
#'   - gtex_op: '`sm config["PROC_RESULTS"] + "/exons_by_gene_op_v19.Rds"`'
#'   - gencode_op: '`sm config["PROC_RESULTS"] + "/exons_by_gene_op_v29.Rds"`'
#'  type: script
#'---

#TODO: change name of output files for later input! rds vs. Rds

saveRDS(snakemake,  "tmp/count_ranges.snakemake")
suppressPackageStartupMessages({
    library(GenomicFeatures)
    library(GenomicRanges)
})

gtex_txdb <- loadDb(snakemake@input$gtex_txdb)
gencode_txdb <- loadDb(snakemake@input$gencode_txdb)

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


