### 20210604 klutz

### basic annotations (start, end, none, both) for full fds
createFDSAnnotations <- function(fds, txdb){
  print("loading introns")
  #seqlevelsStyle(fds) <- seqlevelsStyle(txdb)[1]
  introns <- unique(unlist(intronsByTranscript(txdb)))
  # reduce the introns to only the actually expressed introns
  fds_known <- fds[unique(to(findOverlaps(introns, rowRanges(fds, type = "j"), type = "equal"))),]
  anno_introns <- as.data.table(rowRanges(fds_known, type="j"))
  anno_introns <- anno_introns[,.(seqnames, start, end, strand)]
  
  #calculate extra columns with mean/median intron expression count
  #add the new columns
  print("adding median count to introns")
  sampleCounts <- K(fds_known, type = "j")
  anno_introns[, meanCount := rowMeans(sampleCounts)]
  anno_introns[, medianCount := rowMedians(as.matrix(sampleCounts))]
  
  anno_introns_ranges <- makeGRangesFromDataFrame(anno_introns, keep.extra.columns = TRUE)
  
  ### get all fds junctions
  fds_junctions <- rowRanges(fds, type = "j")
  
  ### Do the annotation just for the most used intron (highest median expression)
  print("prepare the overlaps")
  overlap <- findOverlaps(fds_junctions, anno_introns_ranges)
  both <- findOverlaps(fds_junctions, anno_introns_ranges, type = "equal")
  start <- findOverlaps(fds_junctions, anno_introns_ranges, type = "start")
  end <- findOverlaps(fds_junctions, anno_introns_ranges, type = "end")
  
  print("start calculating annotations")
  annotations <- vapply(c(1:length(fds_junctions)), function(i){
    
    if(!(i %in% from(overlap))) return("none") # junction has no overlap
    
    ## over: indices of all introns that overlap with the current junction
    over <- to(overlap)[which(from(overlap) == i)]
    ## idx of the overlapping intron that is most freq expressed (highest median count)
    maxExpr <- over[which.max(anno_introns_ranges[over, ]$medianCount)]
    
    b <- to(both)[which(from(both) == i)]
    if(length(b) > 0 & maxExpr %in% b) return("both")
    
    s <- to(start)[which(from(start) == i)]
    if(length(s) > 0 & maxExpr %in% s) return("start")
    
    e <- to(end)[which(from(end) == i)]
    if(length(e) > 0 & maxExpr %in% e) return("end")
    
    return("none") #overlaps but no start/end match
  }, FUN.VALUE = 'a')
  
  rowRanges(fds)$annotatedJunction <- annotations
  print("annotations done")
  return(fds)
}



