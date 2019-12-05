#'---
#' title: Pipeline
#' author:
#' wb:
#'  params:
#'    - tmpdir: '`sm drop.getTmpDir()`'
#'  input:
#'    - indexFile: '`sm config["htmlOutputPath"] + "/indexNames.txt"`' 
#' output:
#'   html_document:
#'    code_folding: hide
#'    code_download: TRUE
#'---

#+ echo=F
saveRDS(snakemake, file.path(snakemake@params$tmpdir, "Pipeline.snakemake"))
# snakemake <- readRDS(".drop/tmp/Pipeline.snakemake")

#+ echo=FALSE
indexNames <- scan(snakemake@input$indexFile, character(), quote = "", sep='\n')
pipeline_names <- gsub("_index.html", "", indexNames)
pipeline_names <- gsub("-", " ",pipeline_names)
pipeline_names <- toupper(pipeline_names)

links <- paste0('* [', pipeline_names ,"](./", indexNames, '){target="_blank"}')
links <- paste(links, collapse = '\n')


#' 
#' `r links`

