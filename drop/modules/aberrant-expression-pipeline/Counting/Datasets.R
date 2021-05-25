#'---
#' title: Counts Overview
#' author:  mumichae, salazar
#' wb:
#'  log:
#'   - snakemake: '`sm str(tmp_dir / "AE" / "Count_Overview.Rds")`'
#'  input: 
#'   - summaries: '`sm expand(config["htmlOutputPath"] + 
#'                "/AberrantExpression/Counting/{annotation}/Summary_{dataset}.html",
#'                annotation=cfg.genome.getGeneVersions(), dataset=cfg.AE.groups)`'
#' output:
#'   html_document:
#'    code_folding: hide
#'    code_download: TRUE
#'---

saveRDS(snakemake, snakemake@log$snakemake)

# Obtain the annotations and datasets
gene_annotation_names <- names(snakemake@config$geneAnnotation)
datasets <- snakemake@config$aberrantExpression$groups

#+ echo=FALSE, results="asis"
devNull <- sapply(datasets, function(name){
  sapply(gene_annotation_names, function(version){
    cat(paste0(
      "<h1>Dataset: ", name, ", annotation: ", version, "</h1>",
      "<p>",
      "</br>", "<a href='AberrantExpression/Counting/", version, "/Summary_", name, ".html'   >Count Summary</a>",
      "</br>", "</p>"
    ))
  })
})
