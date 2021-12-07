#'---
#' title: VCF-BAM Matching Analysis over All Datasets
#' author: 
#' wb:
#'  log:
#'   - snakemake: '`sm str(tmp_dir / "MAE" / "QC_overview.Rds")`'
#'  input:
#'   - html: '`sm expand(
#'             config["htmlOutputPath"] + "/MonoallelicExpression/QC/{dataset}.html",
#'             dataset=cfg.MAE.qcGroups
#'      )`'
#' output:
#'  html_document
#'---

saveRDS(snakemake, snakemake@log$snakemake)

# Obtain the datasets
datasets <- snakemake@config$mae$qcGroups 

#+ echo=FALSE, results="asis"
devNull <- sapply(datasets, function(name){
    cat(paste0(
      "<h1>Dataset: ", name, "</h1>",
      "<p>",
      "</br>", "<a href='MonoallelicExpression/QC/", name, ".html'   >QC overview</a>",
      "</br>", "</p>"
    ))
})
