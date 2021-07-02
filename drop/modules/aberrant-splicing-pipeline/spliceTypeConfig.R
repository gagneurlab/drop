##--------------------------------------------
## required packages for the aberrantSpliceTypes
message("Load aberrant splice type packages")
suppressPackageStartupMessages({
  library(AnnotationDbi)
  library(rtracklayer) #to import that blacklist file
  library(GenomicRanges)
})


addUTRLabels <- function(junctions_dt, txdb){
  psi_positions <- which(junctions_dt$type != "theta")
  colnames(junctions_dt)[which(names(junctions_dt) == "STRAND")] <- "strand2"
  junctions_gr <- makeGRangesFromDataFrame(junctions_dt[psi_positions])
  seqlevelsStyle(txdb) <- seqlevelsStyle(junctions_gr)
  ### UTR labels based on txdb file
  ### add 5' 3' UTR labels
  print("find UTR overlap")
  threes <- unique(from(findOverlaps(junctions_gr, threeUTRsByTranscript(txdb, use.names = TRUE))))
  fives <- unique(from(findOverlaps(junctions_gr, fiveUTRsByTranscript(txdb, use.names = TRUE))))
  junctions_dt[psi_positions, UTR := "no"]
  #print("UTRSSSSSSSSSS:")
  #print(threes)
  #print(fives)
  junctions_dt[psi_positions[threes], UTR := "3"]
  junctions_dt[psi_positions[fives], UTR := "5"]
  colnames(junctions_dt)[which(names(junctions_dt) == "strand2")] <- "STRAND"
  print("UTR labels done")
  #print(junctions_dt)
  return(junctions_dt)
}


#addBlacklistLabels <- function(junctions_dt, blacklist_gr){
#  ### add the blacklist information
#  print("Set up blacklist file")
#  psi_positions <- which(junctions_dt$type != "theta")
#  colnames(junctions_dt)[which(names(junctions_dt) == "STRAND")] <- "strand2"
#  junctions_gr <- makeGRangesFromDataFrame(junctions_dt[psi_positions])
#  
#  #blacklist_gr =  import(blacklist_file, format = "BED")
#  seqlevelsStyle(blacklist_gr) <- seqlevelsStyle(junctions_gr)
#  
#  ## create overlap with blacklist and annotate extra column
#  print("find blacklist overlap")
#  black_hits <- unique(from(findOverlaps(junctions_gr, blacklist_gr)))
#  junctions_dt[, blacklist := FALSE]
#  junctions_dt[psi_positions[black_hits], blacklist := TRUE]
#  print("blacklist labels done")
#  
#  colnames(junctions_dt)[which(names(junctions_dt) == "strand2")] <- "STRAND"
#  return(junctions_dt)
#}