#'---
#' title: RNA Variant Calling
#' author:
#' wb:
#'  log:
#'    - snakemake: '`sm str(tmp_dir / "RVC" / "Overview.Rds")`'
#'  input:
#'    - done_files : '`sm expand(os.path.join(str(cfg.processedDataDir) + "/rnaVariantCalling/{dataset}_alt{minAlt}_done.txt"),
#'                       dataset = cfg.RVC.groups,minAlt = getMinAlt())`'
#' output:
#'   html_document:
#'    code_folding: hide
#'    code_download: TRUE
#'---

print("RVC Done")
