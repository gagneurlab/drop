#'---
#' title: Genetic Diagnosis - Pipeline steps in Demo Project
#' author: salazar
#' wb:
#'  params:
#'    - tmpdir: '`sm drop.getTmpDir()`'
#'  input:
#'    - indexFile: '`sm parser.getProcDataDir() + "/indexNames.txt"`' 
#' output:
#'   html_document:
#'    code_folding: hide
#'    code_download: TRUE
#'---

#+ echo=F
saveRDS(snakemake, file.path(snakemake@params$tmpdir, "Pipeline.snakemake"))
# snakemake <- readRDS(".drop/tmp/Pipeline.snakemake")

#' # Links for Pipeline Steps
indexNames <- scan(snakemake@input$indexFile, character(), quote = "", sep='\n')
pipeline_names <- gsub("_index.html", "", indexNames)
pipeline_names <- gsub("-", " ",pipeline_names)
pipeline_names <- toupper(pipeline_names)

links <- paste('[', pipeline_names ,"](./", indexNames, '){target="_blank"}', sep = '')
links <- paste(links, sep = '\n')


#' `r links`

