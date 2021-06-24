##--------------------------------------------
## required packages for the aberrantSpliceTypes
message("Load aberrant splice type packages")
suppressPackageStartupMessages({
  library(AnnotationDbi)
  library(rtracklayer) #to import that blacklist file
})