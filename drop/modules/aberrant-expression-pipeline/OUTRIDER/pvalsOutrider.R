#'---
#' title: P value calculation for OUTRIDER
#' author: Ines Scheller
#' wb:
#'  log:
#'   - snakemake: '`sm str(tmp_dir / "AE" / "{annotation}" / "{dataset}" / "pvalsOUTRIDER.Rds")`'
#'  params:
#'   - ids: '`sm lambda w: sa.getIDsByGroup(w.dataset, assay="RNA")`'
#'   - parse_subsets_for_FDR: '`sm str(projectDir / ".drop" / "helpers" / "parse_subsets_for_FDR.R")`'
#'   - genes_to_test: '`sm cfg.AE.get("genesToTest")`'
#'  input:
#'   - ods_fitted: '`sm cfg.getProcessedResultsDir() + 
#'           "/aberrant_expression/{annotation}/outrider/{dataset}/ods_fitted.Rds"`'
#'   - gene_name_mapping: '`sm cfg.getProcessedDataDir() + "/preprocess/{annotation}/gene_name_mapping_{annotation}.tsv"`'
#'  output:
#'   - ods_with_pvals: '`sm cfg.getProcessedResultsDir() + 
#'           "/aberrant_expression/{annotation}/outrider/{dataset}/ods.Rds"`'
#'  type: script
#'  threads: 30
#'---


#+ echo=F
saveRDS(snakemake, snakemake@log$snakemake)

suppressPackageStartupMessages({
    library(OUTRIDER)
    library(SummarizedExperiment)
    library(ggplot2)
    library(data.table)
    library(dplyr)
    library(magrittr)
    library(tools)
    library(yaml)
})

ods <- readRDS(snakemake@input$ods_fitted)
implementation <- snakemake@config$aberrantExpression$implementation
register(MulticoreParam(snakemake@threads))

# read in gene subsets from sample anno if present (returns NULL if not present)
source(snakemake@params$parse_subsets_for_FDR)
outrider_sample_ids <- snakemake@params$ids
subsets <- parse_subsets_for_FDR(snakemake@params$genes_to_test,
                                 sampleIDs=outrider_sample_ids)
subsets_final <- convert_to_geneIDs(subsets, snakemake@input$gene_name_mapping)

# P value calculation
message(date(), ": P-value calculation ...")
ods <- computePvalues(ods, subsets=subsets_final)
message(date(), ": Zscore calculation ...")
ods <- computeZscores(ods, 
                      peerResiduals=grepl('^peer$', implementation))

# save ods with pvalues
saveRDS(ods, snakemake@output$ods_with_pvals)
