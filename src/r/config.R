# R Script
# author: baderda
##############################################################################

# allow group users to read AND write 
Sys.umask(mode = "0002")

##--------------------------------------------
## required packages

suppressPackageStartupMessages({
    library(Biobase)
    library(data.table)
    library(DESeq2)
    library(GenomicFeatures)
    library(knitr)
    library(plotly)
    library(rmarkdown)
    library(Rsamtools)
    library(tidyr)
})

##--------------------------------------------
## Knitr config
opts_knit$set(root.dir = getwd())
#opts_chunk$set(echo=FALSE, cache=F)


##--------------------------------------------
## functions

stopifnot(dir.exists("../gagneurlab_shared"))

source("src/r/functions/load_rscripts_from_folder.R")
load_rscripts_from_folder("src/r/functions/")
load_rscripts_from_folder("../sample_annotation/src/r/functions/")



##--------------------------------------------
## parameters

# folders
DATADIR     <- "/s/project/mitoMultiOmics/"		# main project folder on ouga
RAWDIR      <- file.path(DATADIR,"raw_data/")	# all original files from the collaborators (read only files) 

# results
PROC_RESULTS <- "/s/project/genetic_diagnosis/processed_results/"
PROC_DATA <- "/s/project/genetic_diagnosis/processed_data/"


# files
FILE_GO_HUMAN    <- "/s/genomes/human/GO/gene_association.goa_human"


##--------------------------------------------
# cutoffs

#' ## low expression filter
#'
LOW_EXPR_QUANTILE= 0.95
RELIABLE_PROT_FRACTION_NA = 0.5

PADJ_METHOD <- 'hochberg'
PADJ_LIMIT <- 0.05
ZSCORE_LIMIT <- 3


##--------------------------------------------
## data

SAMPLE_ANNOTATION <- load_sample_annotation()
