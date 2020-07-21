#'---
#' title: Counts Overview
#' author:  mumichae, salazar
#' wb:
#'  params:
#'   - ids: '`sm parser.outrider_ids`'
#'   - tmpdir: '`sm drop.getMethodPath(METHOD, "tmp_dir")`'
#'  input: 
#'   - summaries: '`sm expand(config["htmlOutputPath"] + 
#'                "/AberrantExpression/Counting/{annotation}/Summary_{dataset}.html",
#'                annotation=list(config["geneAnnotation"].keys()), dataset=parser.outrider_ids)`'
#' output:
#'   html_document:
#'    code_folding: hide
#'    code_download: TRUE
#'---

saveRDS(snakemake, file.path(snakemake@params$tmpdir, "counting_overview.snakemake") )
# snakemake <- readRDS(".drop/tmp/AE/counting_overview.snakemake")

# Obtain the annotations and datasets
gene_annotation_names <- names(snakemake@config$geneAnnotation)
datasets <- snakemake@config$aberrantExpression$groups

#+ echo=FALSE, results="asis"
devNull <- sapply(datasets, function(name){
  sapply(gene_annotation_names, function(version){
  cat(paste0(
    "<h1>Dataset: ", name, "</h1>",
    "<p>",
    "</br>", "<a href='AberrantExpression/Counting/", version, "/Summary_", name, ".html'   >Count Summary</a>",
    "</br>", "</p>"
  ))
  })
})
