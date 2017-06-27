#'---
#' title: R script
#' author: Daniel Bader
#' output: 
#'   html_document:
#'     toc_float: yes
#'     toc: yes
#'---
#'



write_tsv <- function(
    robject, file, row.names = FALSE, ...
){
    write.table(
        x=robject, 
        file=file, 
        quote=FALSE, 
        sep='\t', 
        row.names= row.names, ...)
}






