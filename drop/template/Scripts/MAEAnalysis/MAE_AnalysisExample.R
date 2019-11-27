#'---
#' title: Analysis Example
#' author: salazar
#' wb:
#'  input:
#'    - indexFile: '`sm parser.getProcDataDir() + "/indexNames.txt"  `' 
#' output:
#'   html_document:
#'    code_folding: hide
#'    code_download: TRUE
#'---

#+ echo=F
saveRDS(snakemake, paste0(snakemake@config$tmpdir, "/MAE_Analysis.snakemake") )
# snakemake <- readRDS(paste0(snakemake@config$tmpdir, "/MAE_Analysis.snakemake") )

#' ... Put your results, plots and further analysis here :)
 


