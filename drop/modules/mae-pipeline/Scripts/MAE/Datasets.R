#'---
#' title: MAE analysis over all datasets
#' author: 
#' wb:
#'  params:
#'   - tmpdir: '`sm drop.getMethodPath(METHOD, "tmp_dir")`'
#'  input:
#'   - html: '`sm expand(config["htmlOutputPath"] + 
#'             "/MAE/{dataset}--{annotation}_results.html",
#'              annotation=list(config["geneAnnotation"].keys()),
#'              dataset=config["mae"]["groups"])`'
#' output:
#'  html_document
#'---


saveRDS(snakemake, file.path(snakemake@params$tmpdir, "overview.snakemake") )
# snakemake <- readRDS(".drop/tmp/MAE/overview.snakemake")

# Obtain the annotations and datasets
datasets <- snakemake@config$mae$groups 
gene_annotation_names <- names(snakemake@config$geneAnnotation)

#+ echo=FALSE, results="asis"
devNull <- sapply(datasets, function(name){
  sapply(gene_annotation_names, function(version){
    cat(paste0(
      "<h1>Dataset: ", name, "</h1>",
      "<p>",
      "</br>", "<a href='MAE/", name, "--", version, "_results.html'   >MAE results</a>",
      "</br>", "</p>"
    ))
  })
})

