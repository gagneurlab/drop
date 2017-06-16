##--------------------------------------------
## required packages

## data + convenience libraries
library(data.table)
library(reshape2) ## functons like melt etc by Hadley Wickham
library(tidyr)
library(readr)
library(stringr)
library(Matrix)
library(magrittr)
library(broom)
library(futile.logger)

## plotting libraries
library(ggplot2)
## library(gplots)                 #heatmap.2 function
## library(scales)                 #plot(..., col = alpha(color, 0.5)) function for base r plotting

## Parallel processing
library(parallel)
library(BiocParallel)

## ML libraries
library(caret)
library(glmnet)
library(randomForest)

## Genomic libraries
library(Biostrings)
## library(Rsamtools)		# provides an interface to BAM files
## library(GenomicAlignments)	# Representation and manipulation of short genomic alignment
library(GenomicRanges)		# The GenomicRanges package serves as the foundation for representing genomic
library(rtracklayer)            #read in all the formats
library(GenomicFeatures)
# library(VariantAnnotation)  	# handling vcf/bcf files and some other variant annotations

##--------------------------------------------
## Global variables

STD_CHROMOSOMES <- paste0("chr", c(1:22, "X","Y"))

## PHhome_ase <- "~/ase-grant"
PHhome <- "."
PHgshared <- "../gagneurlab_shared"
PHdata="./data"
PHdatao="./data-offline"
##--------------------------------------------
## User defined functions

## IMPORTANT: local functions have to be loaded last, as I could override some function definitions
## get the function_get_helmholtz_file

## general functions
sapply(list.files(path = file.path(PHgshared,"r/avsec_utils"), pattern=".*\\.R$", full.names=TRUE),source,.GlobalEnv)

## Exome functions
sapply(list.files(path = file.path(PHgshared,"r/exome"), pattern=".*\\.R$", full.names=TRUE),source,.GlobalEnv)

## knitr helper functions
sapply(list.files(path = file.path(PHgshared,"r/knitr_helper"), pattern=".*\\.R$", full.names=TRUE),source,.GlobalEnv)

## kallisto R function
sapply(list.files(path = file.path(PHgshared,"bash/isoform_expression_with_kallisto/r/functions/"),
                  pattern=".*\\.R$", full.names=TRUE),source,.GlobalEnv)


## project specific functions
sapply(list.files(path = file.path(PHhome,"./src/r/functions"), pattern=".*\\.R$",recursive=FALSE,
                  full.names = T),source,.GlobalEnv)
##--------------------------------------------
## Knitr config
library(knitr)
library(rmarkdown)
opts_knit$set(root.dir = getwd())
opts_chunk$set(echo=FALSE, cache=F)
