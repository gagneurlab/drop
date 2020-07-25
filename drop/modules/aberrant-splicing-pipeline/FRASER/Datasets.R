#'---
#' title: Full FRASER analysis over all datasets
#' wb:
#'  log:
#'    - snakemake: '`sm str(tmp_dir / "AS" / "FRASER_datasets.Rds")`'
#'  input:
#'   - fraser_summary: '`sm expand(config["htmlOutputPath"] + 
#'                     "/AberrantSplicing/{dataset}_summary.html", dataset=cfg.AS.groups)`'
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
        "</br>", "<a href='AberrantSplicing/", name, "_summary.html'        >FRASER Summary</a>",
        "</br>", "</p>"
    ))
})
