### 20210702 karoline lutz
### blacklist annotation for aberrant splicing and aberrant expression
### expression junctions dont have type, STRAND 

addBlacklistLabels <- function(junctions_dt, blacklist_gr, aberrant){
  ### add the blacklist information
  print("Set up blacklist file")
  if(aberrant == "splicing"){
    print("Set up aberrant splicing granges")
    
    colnames(junctions_dt)[which(names(junctions_dt) == "STRAND")] <- "strand2"
    
    junctions_gr <- makeGRangesFromDataFrame(junctions_dt)
  }else{
    print("Set up aberrant expression granges")
    junctions_gr <- makeGRangesFromDataFrame(junctions_dt)
  }

  seqlevelsStyle(blacklist_gr) <- seqlevelsStyle(junctions_gr)
  
  ## create overlap with blacklist and annotate extra column
  print("find blacklist overlap")
  black_hits <- unique(from(findOverlaps(junctions_gr, blacklist_gr)))
  junctions_dt[, blacklist := FALSE]
  
  if(aberrant == "splicing"){
    junctions_dt[black_hits, blacklist := TRUE]
    colnames(junctions_dt)[which(names(junctions_dt) == "strand2")] <- "STRAND"
  }else{
    junctions_dt[black_hits, blacklist := TRUE]
  }
  
  print("blacklist labels done")
  
  return(junctions_dt)
}