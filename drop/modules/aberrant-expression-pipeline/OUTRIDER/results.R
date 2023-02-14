#'---
#' title: OUTRIDER Results
#' author: mumichae
#' wb:
#'  log:
#'   - snakemake: '`sm str(tmp_dir / "AE" / "{annotation}" / "{dataset}" / "OUTRIDER_results.Rds")`'
#'  params:
#'   - padjCutoff: '`sm cfg.AE.get("padjCutoff")`'
#'   - zScoreCutoff: '`sm cfg.AE.get("zScoreCutoff")`'
#'   - hpoFile: '`sm cfg.get("hpoFile")`'
#'  input:
#'   - add_HPO_cols: '`sm str(projectDir / ".drop" / "helpers" / "add_HPO_cols.R")`'
#'   - ods: '`sm cfg.getProcessedResultsDir() + "/aberrant_expression/{annotation}/outrider/{dataset}/ods.Rds"`'
#'   - gene_name_mapping: '`sm cfg.getProcessedDataDir() + "/preprocess/{annotation}/gene_name_mapping_{annotation}.tsv"`'
#'   - input_params: '`sm cfg.getProcessedDataDir() + "/aberrant_expression/{annotation}/params/results/{dataset}_resultParams.csv"`'
#'  output:
#'   - results: '`sm cfg.getProcessedResultsDir() + "/aberrant_expression/{annotation}/outrider/{dataset}/OUTRIDER_results.tsv"`'
#'   - results_all: '`sm cfg.getProcessedResultsDir() + "/aberrant_expression/{annotation}/outrider/{dataset}/OUTRIDER_results_all.Rds"`'
#'  type: script
#'---

saveRDS(snakemake, snakemake@log$snakemake)
source(snakemake@input$add_HPO_cols)

suppressPackageStartupMessages({
    library(dplyr)
    library(data.table)
    library(ggplot2)
    library(SummarizedExperiment)
    library(OUTRIDER)
})

ods <- readRDS(snakemake@input$ods)
res <- results(ods, padjCutoff = snakemake@params$padjCutoff,
			   zScoreCutoff = snakemake@params$zScoreCutoff, all = TRUE)

# Add fold change
res[, foldChange := round(2^l2fc, 2)]

# Save all the results and significant ones
saveRDS(res, snakemake@output$results_all)

# Subset to significant results
padj_cols <- grep("padjust", colnames(res), value=TRUE)
res <- res[do.call(pmin, c(res[,padj_cols, with=FALSE], list(na.rm = TRUE))) 
                <= snakemake@params$padjCutoff &
            abs(zScore) >= snakemake@params$zScoreCutoff]

gene_annot_dt <- fread(snakemake@input$gene_name_mapping)
if(!is.null(gene_annot_dt$gene_name)){
  if(grepl('ENSG00', res[1,geneID]) & grepl('ENSG00', gene_annot_dt[1,gene_id])){
    res <- merge(res, gene_annot_dt[, .(gene_id, gene_name)],
                 by.x = 'geneID', by.y = 'gene_id', sort = FALSE, all.x = TRUE)
    setnames(res, 'gene_name', 'hgncSymbol')
    res <- cbind(res[, .(hgncSymbol)], res[, - 'hgncSymbol'])
  }
}

# Add HPO terms, requires online connection and for there to be annotated HPO terms
sa <- fread(snakemake@config$sampleAnnotation, 
              colClasses = c(RNA_ID = 'character', DNA_ID = 'character'))
if(!is.null(sa$HPO_TERMS) & nrow(res) > 0){
  if(!all(is.na(sa$HPO_TERMS)) & ! all(sa$HPO_TERMS == '')){
    res <- add_HPO_cols(res, hpo_file = snakemake@params$hpoFile)
  }
}


# Save results
fwrite(res, snakemake@output$results, sep = "\t", quote = F)
