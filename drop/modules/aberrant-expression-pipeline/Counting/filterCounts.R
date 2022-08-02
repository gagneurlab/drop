#'---
#' title: Filter Counts for OUTRIDER
#' author: Michaela Mueller
#' wb:
#'  log:
#'   - snakemake: '`sm str(tmp_dir / "AE" / "{annotation}" / "{dataset}" / "filter.Rds")`'
#'  input:
#'   - counts: '`sm cfg.getProcessedDataDir() +
#'              "/aberrant_expression/{annotation}/outrider/{dataset}/total_counts.Rds"`'
#'   - txdb: '`sm cfg.getProcessedDataDir() +
#'            "/preprocess/{annotation}/txdb.db"`'
#'  output:
#'   - ods: '`sm cfg.getProcessedResultsDir() +
#'           "/aberrant_expression/{annotation}/outrider/{dataset}/ods_unfitted.Rds"`'
#'  type: script
#'---

saveRDS(snakemake, snakemake@log$snakemake)

suppressPackageStartupMessages({
    library(data.table)
    library(GenomicFeatures)
    library(SummarizedExperiment)
    library(OUTRIDER)
})

counts <- readRDS(snakemake@input$counts)
ods <- OutriderDataSet(counts)
txdb <- loadDb(snakemake@input$txdb)

# filter not expressed genes
fpkmCutoff <- snakemake@config$aberrantExpression$fpkmCutoff
ods <- filterExpression(ods, gtfFile=txdb, filter=FALSE,
                        fpkmCutoff=fpkmCutoff, addExpressedGenes=TRUE)

# add column for genes with at least 1 gene
rowData(ods)$counted1sample = rowSums(assay(ods)) > 0

# External data check
if (is.null(ods@colData$GENE_COUNTS_FILE)){ #column does not exist in sample annotation table
    has_external <- FALSE
}else if(all(is.na(ods@colData$GENE_COUNTS_FILE))){ #column exists but it has no values
    has_external <- FALSE
}else if(all(ods@colData$GENE_COUNTS_FILE == "")){ #column exists with non-NA values but this group has all empty strings
    has_external <- FALSE
}else{ #column exists with non-NA values and this group has at least 1 non-empty string
    has_external <- TRUE
}

if(has_external){
    ods@colData$isExternal <- as.factor(ods@colData$GENE_COUNTS_FILE != "")
}else{
    ods@colData$isExternal <- as.factor(FALSE)
}


# Save the ods before filtering to preserve the original number of genes
saveRDS(ods, snakemake@output$ods)
