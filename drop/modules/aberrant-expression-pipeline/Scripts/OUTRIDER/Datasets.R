#'---
#' title: Results Overview
#' author: mumichae
#' wb:
#'  params:
#'    - tmpdir: '`sm drop.getMethodPath(METHOD, "tmp_dir")`'
#'  input:
#'    - summaries: '`sm expand(config["htmlOutputPath"] + 
#'                 "/AberrantExpression/Outrider/{annotation}/Summary_{dataset}.html",
#'                 annotation=config["geneAnnotation"].keys() , dataset=config["aberrantExpression"]["groups"])`'
#' output:
#'   html_document:
#'    code_folding: hide
#'    code_download: TRUE
#'---

saveRDS(snakemake, file.path(snakemake@params$tmpdir, "outrider_overview.snakemake"))
# snakemake <- readRDS(".drop/tmp/AE/outrider_overview.snakemake")

# Obtain the annotations and datasets
datasets <- snakemake@config$aberrantExpression$groups 
gene_annotation_names <- names(snakemake@config$geneAnnotation)

#+ echo=FALSE, results="asis"
devNull <- sapply(datasets, function(name){
  sapply(gene_annotation_names, function(version){
    cat(paste0(
      "<h1>Dataset: ", name, "</h1>",
      "<p>",
      "</br>", "<a href='AberrantExpression/Outrider/", version, "/Summary_", name, ".html'   >OUTRIDER Summary</a>",
      "</br>", "</p>"
    ))
  })
})
