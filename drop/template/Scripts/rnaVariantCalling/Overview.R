#'---
#' title: RNA Variant Calling
#' author: Nick Smith
#' wb:
#'  log:
#'    - snakemake: '`sm str(tmp_dir / "RVC" / "Overview.Rds")`'
#'  input:
#'    - vcfFilesBatch: '`sm expand(cfg.getProcessedDataDir() +
#'                      "/rnaVariantCalling/out/all_samples_haplocaller/" + 
#'                      "{dataset}/{dataset}.processed.vcf.gz",
#'                  dataset=cfg.RVC.groups)`'
#'    - vcfFilesMasked: '`sm expand(cfg.getProcessedDataDir() +
#'                      "/rnaVariantCalling/out/sample_haplocaller/" + 
#'                      "{sample}/{sample}.vcf.gz",
#'                  sample=cfg.RVC.batchIDs)`'
#' output:
#'   html_document:
#'    code_folding: hide
#'    code_download: TRUE
#'---

#+ echo=FALSE
library(data.table)
library(ggplot2)

saveRDS(snakemake, snakemake@log$snakemake)

