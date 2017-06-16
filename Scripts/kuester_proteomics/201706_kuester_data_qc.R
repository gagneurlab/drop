#'---
#' title: R script
#' author: Daniel Bader
#' wb:
#'   input: 
#'   output: 
#' output: 
#'   html_document:
#'     toc_float: yes
#'     toc: yes
#'---
#'

source("src/r/config.R")


protdir <- file.path(RAWDIR, "proteome", "20170614_kopajtich_kuester_proteome")
files_kuester <- list.files(protdir, pattern = "^m.*txt$", full.names = T)

#' Kuester proteome files:
print(files_kuester)


