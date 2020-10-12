#'---
#' title: FRASER counting analysis over all datasets
#' wb:
#'  log:
#'    - snakemake: '`sm str(tmp_dir / "AS" / "CountingOverview.Rds")`'
#'  input:
#'   - counting_summary: '`sm expand(config["htmlOutputPath"] + 
#'                     "/AberrantSplicing/{dataset}_countSummary.html",
#'                     dataset=cfg.AS.groups)`'
#' output:
#'  html_document
#'---

saveRDS(snakemake, snakemake@log$snakemake)

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
