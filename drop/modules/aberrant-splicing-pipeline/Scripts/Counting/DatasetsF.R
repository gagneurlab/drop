#'---
#' title: FRASER counting analysis over all datasets
#' wb:
#'  params:
#'   - tmpdir: '`sm drop.getMethodPath(METHOD, "tmp_dir")`'
#'  input:
#'   - counting_summary: '`sm expand(config["htmlOutputPath"] + 
#'                     "/AberrantSplicing/{dataset}_countSummary.html",
#'                     dataset=config["aberrantSplicing"]["groups"])`'
#' output:
#'  html_document
#'---

saveRDS(snakemake, file.path(snakemake@params$tmpdir, "FRASER_cs.snakemake"))
# snakemake <- readRDS(".drop/tmp/AS/FRASER_cs.snakemake")

datasets <- snakemake@config$aberrantSplicing$groups

#+ echo=FALSE, results="asis"
devNull <- sapply(datasets, function(name){
  cat(paste0(
    "<h1>Dataset: ", name, "</h1>",
    "<p>",
    "</br>", "<a href='AberrantSplicing/", name, "_countSummary.html'   >Count Summary</a>",
    "</br>", "</p>"
  ))
})