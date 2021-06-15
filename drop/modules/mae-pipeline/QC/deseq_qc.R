#'---
#' title: MAE test on qc variants
#' author: vyepez
#' wb:
#'  log:
#'   - snakemake: '`sm str(tmp_dir / "MAE" / "deseq" / "QC--{rna}.Rds")`'
#'  input:
#'   - qc_counts: '`sm cfg.getProcessedDataDir() + "/mae/allelic_counts/QC--{rna}.csv.gz" `'
#'  output:
#'   - mae_res: '`sm cfg.getProcessedDataDir() + "/mae/RNA_GT/{rna}.Rds"`'
#'  type: script
#'---

saveRDS(snakemake, snakemake@log$snakemake)

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(tMAE)
})

# Read MA counts for qc
qc_counts <- fread(snakemake@input$qc_counts, fill=TRUE)
# qc_counts <- qc_counts[!is.na(position)]

# Run DESeq
rmae <- DESeq4MAE(qc_counts, minCoverage = 10)
rmae[, RNA_GT := '0/1']
rmae[altRatio < .2, RNA_GT := '0/0']
rmae[altRatio > .8, RNA_GT := '1/1']
rmae[, position := as.numeric(position)]

# Convert to granges
qc_gr <- GRanges(seqnames = rmae$contig, 
                 ranges = IRanges(start = rmae$position, end = rmae$position), 
                 strand = '*')
mcols(qc_gr) = DataFrame(RNA_GT = rmae$RNA_GT)

saveRDS(qc_gr, snakemake@output$mae_res)
