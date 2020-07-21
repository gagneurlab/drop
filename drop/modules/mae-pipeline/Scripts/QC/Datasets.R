#'---
#' title: VCF-BAM Matching Analysis over All Datasets
#' author: 
#' wb:
#'  params:
#'   - tmpdir: '`sm drop.getMethodPath(METHOD, "tmp_dir")`'
#'  input:
#'   - html: '`sm expand(config["htmlOutputPath"] + "/QC/{dataset}.html",
#'             dataset=config["mae"]["qcGroups"])`'
#' output:
#'  html_document
#'---

saveRDS(snakemake, file.path(snakemake@params$tmpdir, "overview_qc.snakemake") )
# snakemake <- readRDS(".drop/tmp/MAE/overview_qc.snakemake")

# Obtain the datasets
datasets <- snakemake@config$mae$qcGroups 

#+ echo=FALSE, results="asis"
devNull <- sapply(datasets, function(name){
    cat(paste0(
      "<h1>Dataset: ", name, "</h1>",
      "<p>",
      "</br>", "<a href='QC/", name, ".html'   >QC overview</a>",
      "</br>", "</p>"
    ))
})
