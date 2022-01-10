#'---
#' title: RNA Variant Calling
#' author: Nick Smith
#' wb:
#'  log:
#'    - snakemake: '`sm str(tmp_dir / "RVC" / "Overview.Rds")`'
#'  input:
#'    - singleVCF: '`sm createSingleVCF() `'
#'    - annotatedVCF: '`sm expand(os.path.join(
#'                        cfg.processedResultsDir,
#'                        "rnaVariantCalling/out/batch_vcfs", "{dataset}",
#'                        "{dataset}_{annotation}.annotated.vcf.gz"), 
#'                    annotation = cfg.get("geneAnnotation"), dataset = cfg.RVC.groups) `'
#' output:
#'   html_document:
#'    code_folding: hide
#'    code_download: TRUE
#'---

#+ echo=FALSE
library(data.table)
library(ggplot2)

saveRDS(snakemake, snakemake@log$snakemake)

