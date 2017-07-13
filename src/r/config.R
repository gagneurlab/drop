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
    library(knitr)
    library(plotly)
    library(rmarkdown)
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
NGSDIR      <- file.path(RAWDIR, "helmholtz/")	# all raw sequencing data (BAM,VCF,RNAseq,...)
BADERDIR    <- file.path(DATADIR,'tmp_baderda/')

# results
#PROCDIR     <- file.path(DATADIR,'processed_expression/')	# folder for the processed data
TIDYDIR <- "/s/project/patient_report/tidy_results/"


# files
FILE_GO_HUMAN    <- "/s/genomes/human/GO/gene_association.goa_human"


##--------------------------------------------
# cutoffs

#' ## low expression filter
#'
LOW_EXPR_QUANTILE= 0.95
RELIABLE_PROT_FRACTION_NA = 0.5

PADJ_METHOD <- 'hochberg'


##--------------------------------------------
## data

SAMPLE_ANNOTATION <- load_sample_annotation()
