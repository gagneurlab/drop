##--------------------------------------------
## required packages

suppressPackageStartupMessages({
    library(data.table)
    library(knitr)
    library(rmarkdown)
})

##--------------------------------------------
## Knitr config
opts_knit$set(root.dir = getwd())
opts_chunk$set(echo=FALSE, cache=F)


##--------------------------------------------
## functions

stopifnot(dir.exists("../gagneurlab_shared"))

source("src/r/functions/load_rscripts_from_folder.R")
load_rscripts_from_folder("src/r/functions/")



