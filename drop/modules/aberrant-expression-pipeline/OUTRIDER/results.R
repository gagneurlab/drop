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
#'   - gene_name_mapping: '`sm cfg.getProcessedDataDir() + "/aberrant_expression/{annotation}/gene_name_mapping_{annotation}.tsv"`'
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
res <- results(ods, all = TRUE)

# Add fold change
res[, foldChange := round(2^l2fc, 2)]

# Save all the results and significant ones
saveRDS(res, snakemake@output$results_all)

# Subset to significant results
res <- res[padjust <= snakemake@params$padjCutoff &
               abs(zScore) > snakemake@params$zScoreCutoff]

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
sa <- fread(snakemake@config$sampleAnnotation)
if(!is.null(sa$HPO_TERMS) & nrow(res) > 0){
  if(!all(is.na(sa$HPO_TERMS)) & ! all(sa$HPO_TERMS == '')){
    res <- add_HPO_cols(res, hpo_file = snakemake@params$hpoFile)
  }
}


# Save results 
fwrite(res, snakemake@output$results, sep = "\t", quote = F)

web_dir <- snakemake@config$webDir
if (!is.null(web_dir)) {
    pub_res <- paste0(web_dir, 
                      "/aberrant_expression/results/",{snakemake@wildcards$annotation},"/outrider/",
                      {snakemake@wildcards$dataset},"/OUTRIDER_results.tsv")
    fwrite(res, pub_res, sep = "\t", quote = F)
}
