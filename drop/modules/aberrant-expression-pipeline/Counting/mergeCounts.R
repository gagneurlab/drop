#'---
#' title: Merge the counts for all samples
#' author: Michaela Müller
#' wb:
#'  log:
#'   - snakemake: '`sm str(tmp_dir / "AE" / "{annotation}" / "{dataset}" / "merge.Rds")`'
#'  params:
#'   - exCountIDs: '`sm lambda w: sa.getIDsByGroup(w.dataset, assay="GENE_COUNT")`'
#'  input: 
#'   - counts: '`sm lambda w: cfg.AE.getCountFiles(w.annotation, w.dataset)`'
#'   - count_ranges: '`sm cfg.getProcessedDataDir() + "/aberrant_expression/{annotation}/count_ranges.Rds" `'
#'   - input_params: '`sm cfg.getProcessedDataDir() + "/aberrant_expression/{annotation}/params/merge/{dataset}_mergeParams.csv"`'
#'  output:
#'    - counts: '`sm cfg.getProcessedDataDir() +
#'               "/aberrant_expression/{annotation}/outrider/{dataset}/total_counts.Rds"`'
#'  threads: 30
#'  type: script
#'---

saveRDS(snakemake, snakemake@log$snakemake)

suppressPackageStartupMessages({
    library(data.table)
    library(dplyr)
    library(BiocParallel)
    library(SummarizedExperiment)
})

register(MulticoreParam(snakemake@threads))
count_ranges <- readRDS(snakemake@input$count_ranges)

# Read counts
counts_list <- bplapply(snakemake@input$counts, function(f){
    if(grepl('Rds$', f))
        assay(readRDS(f))
    else {
        ex_counts <- as.matrix(fread(f), rownames = "geneID")
        print(head(ex_counts))
        stopifnot(! snakemake@params$exCountIDs %in% names(ex_counts))
        subset(ex_counts, select = snakemake@params$exCountIDs)
    }
})
message(paste("read", length(counts_list), 'files'))

# check rownames and proceed only if they are the same
row_names_objects <- lapply(counts_list, rownames)
if( length(unique(row_names_objects)) > 1 ){
  stop('The rows (genes) of the count matrices to be merged are not the same.')
}

# merge counts
merged_assays <- do.call(cbind, counts_list)
total_counts <- SummarizedExperiment(assays=list(counts=merged_assays))
colnames(total_counts) <- gsub('.bam', '', colnames(total_counts))

# assign ranges
rowRanges(total_counts) <- count_ranges

# Add sample annotation data (colData)
sample_anno <- fread(snakemake@config$sampleAnnotation,
                    colClasses = c(RNA_ID = 'character', DNA_ID = 'character'))
sample_anno <- sample_anno[, .SD[1], by = RNA_ID]
col_data <- data.table(RNA_ID = as.character(colnames(total_counts)))
col_data <- left_join(col_data, sample_anno, by = "RNA_ID")
rownames(col_data) <- col_data$RNA_ID
colData(total_counts) <- as(col_data, "DataFrame")
rownames(colData(total_counts)) <- colData(total_counts)$RNA_ID

# save in RDS format
saveRDS(total_counts, snakemake@output$counts)
