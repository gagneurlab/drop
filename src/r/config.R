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
    library(dplyr)
    library(GenomicFeatures)
    library(heatmaply)
    library(knitr)
    library(magrittr)
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

stopifnot(dir.exists("../sample_annotation"))

source("src/r/functions/load_rscripts_from_folder.R")
load_rscripts_from_folder("src/r/functions/")
load_rscripts_from_folder("../sample_annotation/src/r/functions/")



##--------------------------------------------
## parameters

# folders
DATADIR <- "/s/project/mitoMultiOmics"		# main project folder on ouga
RAWDIR  <- file.path(DATADIR,"raw_data")	# all original files from the collaborators (read only files) 

# result folders
PROC_RESULTS <- "/s/project/genetic_diagnosis/processed_results"
PROC_DATA    <- "/s/project/genetic_diagnosis/processed_data"

# files
FILE_GO_HUMAN <- "/s/genomes/human/GO/gene_association.goa_human"


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

## Needed for get_helmholtz_file()
HHZ_FILE_SOURCES =  c("GENOME", "EXOME", "RNA") 
HHZ_FILES = structure(list(BAM = "merged.bam", BAM.BAI = "merged.bam.bai", 
                           BAM.STAR = "STAR.bam", BAM.BAI.STAR = "STAR.bam.bai", VCF = "ontarget.varfilter.dbSNP.plus.checked.vcf", 
                           VCF = "ontarget.varfilter.dbSNP.plus.checked.vcf.gz", VCF.DD = "ontarget.dd.out.vcf", 
                           VCF.DD2 = "ontarget.dd2.out.vcf", VCF.CSI = "ontarget.varfilter.dbSNP.plus.checked.vcf.gz.csi", 
                           VCF.PIN = "pindel.vcf", GVCF = "all.gatk.ontarget.haplotypecaller.gvcf", 
                           GATKVCF = "gatk.ontarget.haplotypecaller.filtered.dbSNP.plus.checked.vcf", 
                           VEP = "_annotated_data_table.rds"), .Names = c("BAM", "BAM.BAI", 
                                                                          "BAM.STAR", "BAM.BAI.STAR", "VCF", "VCF", "VCF.DD", "VCF.DD2", 
                                                                          "VCF.CSI", "VCF.PIN", "GVCF", "GATKVCF", "VEP"))


