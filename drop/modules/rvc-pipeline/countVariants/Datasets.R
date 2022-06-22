#'---
#' title: RVC datasets
#' author: nickhsmith
#' wb:
#'  log:
#'   - snakemake: '`sm str(tmp_dir / "RVC" / "RVC_Datasets.Rds")`'
#'  input:
#'    - summaries: '`sm expand(config["htmlOutputPath"] + 
#'                 "/rnaVariantCalling/{annotation}/Summary_{dataset}.html",
#'                 annotation=cfg.genome.getGeneVersions(), dataset=cfg.RVC.groups)`'
#' output:
#'   html_document:
#'    code_folding: hide
#'    code_download: TRUE
#'---

saveRDS(snakemake, snakemake@log$snakemake)

# Obtain the annotations and datasets
datasets <- snakemake@config$rnaVariantCalling$groups 
gene_annotation_names <- names(snakemake@config$geneAnnotation)

#+ echo=FALSE, results="asis"
devNull <- sapply(datasets, function(name){
  sapply(gene_annotation_names, function(version){
    cat(paste0(
      "<h4>Dataset: ", name, ", annotation: ", version, "</h4>",
      "<p>",
      "</br>", "<a href='rnaVariantCalling/", version, "/Summary_", name, ".html'   >Summary</a>",
      "</br>", "</p>"
    ))
  })
})
