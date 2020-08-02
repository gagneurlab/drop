#'---
#' title: Export counts in tsv format
#' author: Michaela Mueller, vyepez
#' wb:
#'  log:
#'   - snakemake: '`sm str(tmp_dir / "AE" / "{annotation}" / "{dataset}" / "export_{genomeAssembly}.Rds")`'
#'  input: 
#'    - counts: '`sm cfg.getProcessedDataDir() +
#'               "/aberrant_expression/{annotation}/outrider/{dataset}/total_counts.Rds"`'
#'  output:
#'    - export: '`sm cfg.exportCounts.getFilePattern(str_=False) / "geneCounts.tsv.gz"`'
#'  type: script
#'---

saveRDS(snakemake, snakemake@log$snakemake)

suppressPackageStartupMessages({
    library(data.table)
    library(SummarizedExperiment)
})

total_counts <- readRDS(snakemake@input$counts)

# save in exportable format
fwrite(as.data.table(assay(total_counts), keep.rownames = 'geneID'),
       file = snakemake@output$export,
       quote = FALSE, row.names = FALSE, sep = '\t', compress = 'gzip')
