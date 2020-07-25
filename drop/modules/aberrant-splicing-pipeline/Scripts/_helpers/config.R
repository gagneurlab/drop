##--------------------------------------------
## required packages
message("Load packages")
suppressPackageStartupMessages({
    library(markdown)
    library(knitr)
    library(devtools)
    library(yaml)
    library(BBmisc)
    library(GenomicAlignments)
    library(tidyr)
    library(data.table)
    library(dplyr)
    library(plotly)
    library(DelayedMatrixStats)
    library(FRASER)
  library(rhdf5)
})


## helper functions
write_tsv <- function(x, file, row.names = FALSE, ...){
  write.table(x=x, file=file, quote=FALSE, sep='\t', row.names= row.names, ...)
}

extract_params <- function(params) {
    unlist(params)[1]
}

options("FRASER.maxSamplesNoHDF5"=1)
options("FRASER.maxJunctionsNoHDF5"=-1)

h5disableFileLocking()
