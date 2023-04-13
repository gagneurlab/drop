##--------------------------------------------
## required packages
message("Load packages")
suppressPackageStartupMessages({
    library(rmarkdown)
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

options("FRASER.maxSamplesNoHDF5"=0)
options("FRASER.maxJunctionsNoHDF5"=-1)

h5disableFileLocking()

# set psiTypes to run based on preference in config.yaml
cfg <- yaml::read_yaml("config.yaml")
if(cfg$aberrantSplicing$FRASER_version == "FRASER2"){
    pseudocount(0.1)
    psiTypes <- c("jaccard")
    psiTypesNotUsed <- c("psi5", "psi3", "theta")
} else{
    pseudocount(1)
    psiTypes <- c("psi5", "psi3", "theta")
    psiTypesNotUsed <- c("jaccard")
}
