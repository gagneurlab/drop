#'---
#' title: MAE analysis over all datasets
#' author: 
#' wb:
#'  log:
#'    - snakemake: '`sm str(tmp_dir / "MAE" / "overview.Rds")`'
#'  input:
#'    - html: '`sm expand(config["htmlOutputPath"] + 
#'             "/MonoallelicExpression/{dataset}--{annotation}_results.html",
#'              annotation=cfg.genome.getGeneVersions(), dataset=cfg.MAE.groups)`'
#' output:
#'  html_document
#'---

saveRDS(snakemake, snakemake@log$snakemake)

# Obtain the annotations and datasets
datasets <- snakemake@config$mae$groups 
gene_annotation_names <- names(snakemake@config$geneAnnotation)

#+ echo=FALSE, results="asis"
devNull <- sapply(datasets, function(name){
  sapply(gene_annotation_names, function(version){
    cat(paste0(
      "<h1>Dataset: ", name, ", annotation ", version, "</h1>",
      "<p>",
      "</br>", "<a href='MonoallelicExpression/", name, "--", version, "_results.html'   >MAE results</a>",
      "</br>", "</p>"
    ))
  })
})

