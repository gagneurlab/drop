### 20210604 klutz

### basic annotations (start, end, none, both) for full fds

testFct <- function(fds){
  message("start test function")
  test_vector <- rep("test", times = length(rowRanges(fds, type="j")))
  rowRanges(fds)$test = test_vector
  message("end test functiont")
  return(fds)
}

createFDSAnnotations <- function(fds, txdb){
  print("loading introns")
  #seqlevelsStyle(fds) <- seqlevelsStyle(txdb)[1]
  introns <- unique(unlist(intronsByTranscript(txdb)))
  # reduce the introns to only the actually expressed introns
  fds_known <- fds[unique(to(findOverlaps(introns, rowRanges(fds, type = "j"), type = "equal"))),]
  grAnno <- rowRanges(fds_known, type="j")
  anno_introns <- as.data.table(grAnno)
  anno_introns <- anno_introns[,.(seqnames, start, end, strand)]
  
  #calculate extra columns with mean/median intron expression count
  #add the new columns
  print("adding median count to introns")
  sampleCounts <- K(fds_known, type = "j")
  anno_introns[, "meanCount" := rowMeans(sampleCounts)]
  anno_introns[, "medianCount" := rowMedians(as.matrix(sampleCounts))]
  
  anno_introns_ranges <- makeGRangesFromDataFrame(anno_introns, keep.extra.columns = TRUE)
  
  ### get all fds junctions
  fds_junctions <- rowRanges(fds, type = "j")
  
  ### Do the annotation just for the most used intron (highest median expression)
  print("start calculating annotations")
  annotations <- sapply(c(1:length(fds_junctions)), function(i){
    #print(i)
    #print("-------------")
    overlap <- to(findOverlaps(fds_junctions[i], anno_introns_ranges))
    if(length(overlap) == 0) return("none") #no overlap with any intron
    
    expre <- sapply(overlap, function(j){
      elementMetadata(anno_introns_ranges[j])$medianCount
    })
    maxExpr <- which.max(expre)
    
    hit_equal <- from(findOverlaps(fds_junctions[i], anno_introns_ranges[overlap[maxExpr]], type="equal"))
    if(length(hit_equal) > 0) return("both")
    
    hit_start <- from(findOverlaps(fds_junctions[i], anno_introns_ranges[overlap[maxExpr]], type="start"))
    if(length(hit_start) > 0) return("start")
    hit_end   <- from(findOverlaps(fds_junctions[i], anno_introns_ranges[overlap[maxExpr]], type="end"))
    if(length(hit_end) > 0) return("end")
    
    return("none") #overlaps but no start/end match
  })
  
  #table(annotations)
  rowRanges(fds)$annotatedJunction = annotations
  print("annotations done")
  return(fds)
}













