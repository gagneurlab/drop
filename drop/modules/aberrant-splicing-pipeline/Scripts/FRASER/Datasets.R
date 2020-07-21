#'---
#' title: Full FRASER analysis over all datasets
#' wb:
#'  params:
#'   - tmpdir: '`sm drop.getMethodPath(METHOD, "tmp_dir")`'
#'  input:
#'   - fraser_summary: '`sm expand(config["htmlOutputPath"] + 
#'                     "/AberrantSplicing/{dataset}_summary.html", 
#'                     dataset=config["aberrantSplicing"]["groups"])`'
#' output:
#'  html_document
#'---

saveRDS(snakemake, file.path(snakemake@params$tmpdir, "FRASER_99.snakemake"))
# snakemake <- readRDS(".drop/tmp/AS/FRASER_99.snakemake")

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
