### 20210427 karoline lutz
### calculate splice type for aberrant splicing junctions and frameshift


aberrantSpliceType <- function(junctions_dt, fds, txdb){
  print("preparing..")
  psi_positions <- which(junctions_dt$type != "theta")
  colnames(junctions_dt)[which(names(junctions_dt) == "STRAND")] <- "strand2"
  junctions_gr <- makeGRangesFromDataFrame(junctions_dt[psi_positions], keep.extra.columns = T)
  seqlevelsStyle(txdb) <- seqlevelsStyle(junctions_gr)
  
  introns_tmp <- unique(unlist(intronsByTranscript(txdb)))
  exons <- exons(txdb)
  
  #seqlevelsStyle(fds) <- seqlevelsStyle(txdb)[1]
  fds_known <- fds[unique(to(findOverlaps(introns_tmp, rowRanges(fds, type = "j"), type = "equal"))),]
  grIntrons <- rowRanges(fds_known, type="j")
  introns <- as.data.table(grIntrons)
  introns <- introns[,.(seqnames, start, end, strand)]
  
  sampleCounts <- K(fds_known, type = "j")
  introns[, "meanCount" := rowMeans(sampleCounts)]
  introns[, "medianCount" := rowMedians(as.matrix(sampleCounts))]
  intron_ranges <- makeGRangesFromDataFrame(introns, keep.extra.columns = TRUE)
  
  ## prepare the results column
  junctions_dt[, aberrantSpliceType := "NA"]
  junctions_dt[, causesFrameshift := "NA"]
  junctions_dt[which(junctions_dt$annotatedJunction=="both"), aberrantSpliceType := "sameIntron"]
  junctions_dt[which(junctions_dt$annotatedJunction=="both"), causesFrameshift := "unlikely"]
  
  starts <- which(junctions_dt[psi_positions]$annotatedJunction=="start")
  ends <- which(junctions_dt[psi_positions]$annotatedJunction=="end")
  nones <- which(junctions_dt[psi_positions]$annotatedJunction=="none")
  
  print("calculating aberrant splice types")
  print("start junctions")
  start_results <- sapply(starts, function(i){
    #print("---------")
    #print(i)
    ## find the most freq intron that overlaps again
    overlap <- to(findOverlaps(junctions_gr[i], intron_ranges, type = "start"))
    #print(overlap)
    expre <- sapply(overlap, function(j){
      elementMetadata(intron_ranges[j])$medianCount
    })
    maxExpr <- which.max(expre)
    #print(maxExpr)
    
    return(compareEnds(junctions_gr, i, overlap[maxExpr], F, intron_ranges, exons))
  })
  junctions_dt[psi_positions[starts], causesFrameshift:=start_results[2,]]
  junctions_dt[psi_positions[starts], aberrantSpliceType := start_results[1,]]
  
  print("end junctions")
  end_results <- sapply(ends, function(i){
    #print("----------------")
    #print(i)
    ## find the most freq intron that overlaps again
    overlap <- to(findOverlaps(junctions_gr[i], intron_ranges, type = "end"))
    #print(overlap)
    expre <- sapply(overlap, function(j){
      elementMetadata(intron_ranges[j])$medianCount
    })
    maxExpr <- which.max(expre)
    #print(maxExpr)
    
    return(compareStarts(junctions_gr, i, overlap[maxExpr], F, intron_ranges, exons))
    
  })
  junctions_dt[psi_positions[ends], causesFrameshift:=end_results[2,]]
  junctions_dt[psi_positions[ends], aberrantSpliceType := end_results[1,]]
  
  print("none junctions pt1")
  none_results <- sapply(nones, function(i){
    ## find most freq intron
    ## check start and end
    
    #print("---------")
    #print(i)
    ## find the most freq intron that overlaps again
    overlap <- to(findOverlaps(junctions_gr[i], intron_ranges))
    if(length(overlap) == 0) return(c("noOverlap", "inconclusive"))
    
    #print(overlap)
    expre <- sapply(overlap, function(j){
      elementMetadata(intron_ranges[j])$medianCount
    })
    maxExpr <- which.max(expre)
    #print(maxExpr)
    
    ## returns type of exon splicing, frameshift T/F, amount of shift
    st = compareStarts(junctions_gr, i, overlap[maxExpr], T, intron_ranges, exons)
    en = compareEnds(junctions_gr, i, overlap[maxExpr], T, intron_ranges, exons)
    
    ## merge, start and end results
    ## merge exon elongation/truncation
    ## if both likely/unlikely fine
    ## if one is likely -> return likely
    ## if one is notYet -> return notYet
    if((st[1] == "singleExonSkipping" & !(en[1] %in% c("singleExonSkipping", "exonSkipping"))) ||
       (en[1] == "singleExonSkipping" & !(st[1] %in% c("singleExonSkipping", "exonSkipping")))){
      ## only one is single exonSkipping, the other is trunc/elong
      if((as.integer(st[3])+as.integer(en[3]))%%3 != 0){
        frs = "likely"
      }else{ frs = "unlikely"}
      return(c("singleExonSkipping", frs))
    } 
    
    if(st[1] %in% c("exonSkipping", "singleExonSkipping") || en[1] %in% c("exonSkipping", "singleExonSkipping")) return(c("exonSkipping", "inconclusive"))
    
    #print(st)
    #print(en)
    if((as.integer(st[3])+as.integer(en[3]))%%3 != 0){
      frs = "likely"
    }else{ frs = "unlikely"}
    if( st[1] != en[1]){
      #combined = paste(st[1], en[1], sep=',')
      combined = "trunc, elong"
    }else{combined = st[1]}
    #combined = toString(st[1] ) + "," + toString(en[1])
    return(c(combined,frs))
    
  })
  junctions_dt[psi_positions[nones], causesFrameshift:=none_results[2,]]
  junctions_dt[psi_positions[nones], aberrantSpliceType := none_results[1,]]
  
  noLaps <-which(junctions_dt[psi_positions]$aberrantSpliceType=="noOverlap")
  refseq.genes<- genes(txdb)
  
  print("none junctions pt2")
  noLaps_results <- sapply(noLaps, function(i){
    #print("--------------")
    #print(i)
    overlap <- to(findOverlaps(junctions_gr[i], exons))
    ## no overlap with an intron or an exon
    if(length(overlap) == 0) return(checkIntergenic(junctions_gr, i, refseq.genes)) #return(c("noOverlap", "noOverlap"))
    
    ## for the exons, check if splice site is contained in the exon
    for(j in overlap){
      exon_start = start(exons[j])
      exon_end = end(exons[j])
      if(exon_start <= start(junctions_gr[i]) & exon_end >= end(junctions_gr[i])){
        if((end(junctions_gr[i]) - start(junctions_gr[i]) + 1)%%3 != 0){
          frs = "likely"
        }else{ frs = "unlikely"}
        return(c("exonTruncation", frs))
      }
    }
    
    #print(overlap)
    return(c("inconclusive","inconclusive"))
  })
  junctions_dt[psi_positions[noLaps], causesFrameshift:=noLaps_results[2,]]
  junctions_dt[psi_positions[noLaps], aberrantSpliceType := noLaps_results[1,]]
  
  
  ### add distance to closest neighbour gene for intergenic results
  print("adding distances to nearest gene")
  up <- which(junctions_dt[psi_positions]$aberrantSpliceType == "upstream")
  down <- which(junctions_dt[psi_positions]$aberrantSpliceType == "downstream")
  print("Calculate distances")
  
  if(length(up) > 0){
    distanceNearestGene_up <- sapply(up, function(i) min(distance(junctions_gr[i], refseq.genes), na.rm = T))
    if(length(distanceNearestGene_up > 0)){
      junctions_dt[psi_positions[up], distNearestGene := distanceNearestGene_up]
    } else{
      junctions_dt[psi_positions[up], distNearestGene := NA]
      print("No distances found for upstream")
    }
  }else{print("No upstream targets")}
  
  if(length(down) > 0){
    distanceNearestGene_down <- sapply(down, function(i) min(distance(junctions_gr[i], refseq.genes), na.rm = T))
    if(length(distanceNearestGene_down > 0)){
      junctions_dt[psi_positions[down], distNearestGene := distanceNearestGene_down]
    }else{
      junctions_dt[psi_positions[down], distNearestGene := NA]
      print("No distances found for downstream")
    }
  }else{print("No downstream targets")}
  
  colnames(junctions_dt)[which(names(junctions_dt) == "strand2")] <- "STRAND"
  print("done calculating aberrant splice types")
  return(junctions_dt)
}




############# FUNCTIONS ##############################

compareStarts <- function(junctions_gr, i, max_lap, shift_needed, intron_ranges, exons){
  intron_start = start(intron_ranges[max_lap])
  ss_start = start(junctions_gr[i])
  
  ## found the most freq intron with same end again
  ## check if intron starts before splice site -> exon elongation -> FRS -> done
  if(intron_start < ss_start){
    #print(intron_ranges[overlap[maxExpr]])
    #print(test_set[i])
    if(((ss_start - intron_start)%%3) != 0){
      frs = "likely"
    }else{ frs = "unlikely"}
    
    ifelse(shift_needed, return(c("exonElongation", frs, (ss_start - intron_start))), return(c("exonElongation", frs)))
    #return(c("exon elongation", frs))
  }
  
  ## check if splice site ends in following exon -> exon truncation -> FRS -> done
  if(intron_start > ss_start){
    #print(intron_ranges[overlap[maxExpr]])
    #print(test_set[i])
    
    ## create dummy exon find all exons starting from that intron end
    dummy_exon <- GRanges(
      seqnames = toString(seqnames(intron_ranges[max_lap])),
      ranges = IRanges(intron_start-2, end = intron_start -1),
      strand = toString(strand(intron_ranges[max_lap]))
    )
    exonChoices <- to(findOverlaps(dummy_exon, exons, type = "end"))
    #print(exonChoices)
    for(j in exonChoices){
      exon_start = start(exons[j])
      #print(exon_end)
      if(exon_start < ss_start){
        if((end(exons[j]) - ss_start + 1)%%3 != 0){
          frs = "likely"
        }else{frs = "unlikely"}
        ifelse(shift_needed, return(c("exonTruncation", frs, (end(exons[j]) - ss_start + 1))), return(c("exonTruncation", frs)))
        #return(c("exon truncation",frs))
      }
    }
    
    
    
    ## check for single exon skipping
    if(length(exonChoices) == 1){
      
      ## check if there is no other exon within the first intron: splice site end until exon end
      dummyFirstItr <- GRanges(
        seqnames = toString(seqnames(intron_ranges[max_lap])),
        ranges = IRanges(end(exons[exonChoices[1]]) + 1, end(junctions_gr[i])),
        strand = toString(strand(intron_ranges[max_lap]))
      )
      
      if(length(findOverlaps(exons, dummyFirstItr, type = "within")) > 0){
        ## another exon is contained within the most freq used intron
        ifelse(shift_needed, return(c("exonSkipping", "inconclusive", 0)), return(c("exonSkipping", "inconclusive")))
      }
      
      
      secondItr <- GRanges(
        seqnames = toString(intron_ranges[max_lap]@seqnames@values),
        strand = toString(intron_ranges[max_lap]@strand@values),
        ranges = IRanges(ss_start, start(exons[exonChoices[1]]) - 1) # end of exon + 1, end of aberrant junction
      )
      secItrChoices <- to(findOverlaps(secondItr, intron_ranges, type = "end"))
      ### only look at most used one
      expre <- sapply(secItrChoices, function(j){
        elementMetadata(intron_ranges[j])$medianCount
      })
      maxExpr <- which.max(expre)
      
      if(length(secItrChoices) == 0){
        ifelse(shift_needed, return(c("exonSkipping", "inconclusive", 0)), return(c("exonSkipping", "inconclusive")))
      }
      
      if(ss_start >= start(intron_ranges[secItrChoices[maxExpr]])){
        ## check if there is no other exon in that range
        if(length(findOverlaps(exons, intron_ranges[secItrChoices[maxExpr]], type = "within")) == 0){
          ## clear exon skipping, only exon is skipped 
          ## calculate frameshift, skipped exon plus possible exon elongation 
          shift = end(exons[exonChoices[1]]) - start(exons[exonChoices[1]]) + 1 + ss_start - start(intron_ranges[secItrChoices[maxExpr]])
          frs = ifelse(shift%%3 == 0,"unlikely","likely")
          ifelse(shift_needed, return(c("singleExonSkipping", "inconclusive", shift)), return(c("singleExonSkipping", frs)))
        }
      }
    } ## single exon skipping end
    
  }
  
  ## splice site longer than one intron + exon -> not defined for now
  ifelse(shift_needed, return(c("exonSkipping", "inconclusive", 0)), return(c("exonSkipping", "inconclusive")))
  #return(c("notYet", "notYet"))
}


compareEnds <- function(junctions_gr, i, max_lap, shift_needed, intron_ranges, exons){
  intron_end = end(intron_ranges[max_lap])
  ss_end = end(junctions_gr[i])
  
  ## found the most freq intron with same start again
  ## check if intron ends after splice site -> exon elongation -> FRS -> done
  if(intron_end > ss_end){
    #print(intron_ranges[overlap[maxExpr]])
    #print(test_set[i])
    if(((intron_end - ss_end)%%3) != 0){
      frs = "likely"
    }else{ frs = "unlikely"}
    
    ifelse(shift_needed, return(c("exonElongation", frs, (intron_end - ss_end))), return(c("exonElongation", frs)))
    #return(c("exon elongation", frs))
  }
  
  ## check if splice site ends in following exon -> exon truncation -> FRS -> done
  if(intron_end < ss_end){
    #print(intron_ranges[overlap[maxExpr]])
    #print(test_set[i])
    
    ## create dummy exon find all exons starting from that intron end
    dummy_exon <- GRanges(
      seqnames = toString(intron_ranges[max_lap]@seqnames@values),
      ranges = IRanges(intron_end + 1, end = intron_end + 2),
      strand = toString(intron_ranges[max_lap]@strand@values)
    )
    exonChoices <- to(findOverlaps(dummy_exon, exons, type = "start"))
    #print(exonChoices)
    for(j in exonChoices){
      exon_end = end(exons[j])
      #print(exon_end)
      if(exon_end > ss_end){
        if((ss_end - start(exons[j]) + 1)%%3 != 0){
          frs = "likely"
        }else{frs = "unlikely"}
        ifelse(shift_needed, return(c("exonTruncation",frs, (ss_end - start(exons[j]) + 1))), return(c("exonTruncation",frs)))
        #return(c("exon truncation",frs))
      }
    }
    
    ## check for single exon skipping
    if(length(exonChoices) == 1){
      
      ## check if there is no other exon within the first intron: splice site end until exon end
      dummyFirstItr <- GRanges(
        seqnames = toString(seqnames(intron_ranges[max_lap])),
        ranges = IRanges(start(junctions_gr[i]), start(exons[exonChoices[1]]) - 1),
        strand = toString(strand(intron_ranges[max_lap]))
      )
      
      if(length(findOverlaps(exons, dummyFirstItr, type = "within")) > 0){
        ## another exon is contained within the most freq used intron
        ifelse(shift_needed, return(c("exonSkipping", "inconclusive", 0)), return(c("exonSkipping", "inconclusive")))
      }
      
      
      secondItr <- GRanges(
        seqnames = toString(intron_ranges[max_lap]@seqnames@values),
        strand = toString(intron_ranges[max_lap]@strand@values),
        ranges = IRanges(end(exons[exonChoices[1]]) + 1, ss_end) # end of exon + 1, end of aberrant junction
      )
      secItrChoices <- to(findOverlaps(secondItr, intron_ranges, type = "start"))
      ### only look at most used one
      expre <- sapply(secItrChoices, function(j){
        elementMetadata(intron_ranges[j])$medianCount
      })
      maxExpr <- which.max(expre)
      
      if(length(secItrChoices) == 0){
        ifelse(shift_needed, return(c("exonSkipping", "inconclusive", 0)), return(c("exonSkipping", "inconclusive")))
      }
      
      if(ss_end <= end(intron_ranges[secItrChoices[maxExpr]])){
        ## check if there is no other exon in that range
        if(length(findOverlaps(exons, intron_ranges[secItrChoices[maxExpr]], type = "within")) == 0){
          ## clear exon skipping, only exon is skipped 
          ## calculate frameshift, skipped exon plus possible exon elongation at end
          shift = end(exons[exonChoices[1]]) - start(exons[exonChoices[1]]) + 1 + end(intron_ranges[secItrChoices[maxExpr]]) - ss_end
          frs = ifelse(shift%%3 == 0,"unlikely","likely")
          ifelse(shift_needed, return(c("singleExonSkipping", "inconclusive", shift)), return(c("singleExonSkipping", frs)))
        }
      }
    } ## single exon skipping end
    
    
  }
  
  ## splice site longer than one intron + exon -> not defined for now
  ifelse(shift_needed, return(c("exonSkipping", "inconclusive", 0)), return(c("exonSkipping", "inconclusive")))
  #return(c("notYet", "notYet"))
}


checkIntergenic <- function(junctions_gr, i, refseq.genes){
  ## check if start > 1000
  ## start - 1000, end + 1000
  start = start(junctions_gr[i]) 
  #ifelse(start > 1000, start = start - 1000, start = 1)
  #if(start > 1000){
  #  start = start - 1000
  #}else{start = 1}
  
  end = end(junctions_gr[i])  #+ 1000
  if(start + 2 < end){
    start = start + 1
    end = end - 1
  }
  
  test_junction <- GRanges(
    seqnames = seqnames(junctions_gr[i]),
    ranges = IRanges(start, end),
    strand = strand(junctions_gr[i])
  )
  
  ## overlap with introns and exon
  ## IGNORE STRANDS? -> decided its not necessary
  #intron_overlap_length = length(to(findOverlaps(test_junction, intron_ranges)))
  #exons_overlap_length = length(to(findOverlaps(test_junction, exons)))
  #if(intron_overlap_length + exons_overlap_length == 0){
  
  
  ### check if distance to nearest is > 1000 -> intergenic
  ### otherwise up/downstream
  dist = min(distance(test_junction, refseq.genes), na.rm = T)
  #elementMetadata(distanceToNearest(test_junction, refseq.genes))$distance
  if(dist > 0){
    #if(dist > 1000){
    #  print("intergenic")
    #  return(c("intergenic", "unlikely"))
    #}else{
    ## find nearest and compare starts
    if(start(refseq.genes[nearest(junctions_gr[i], refseq.genes)]) > start){
      ifelse(strand(junctions_gr[i]) == "+", return(c("upstream", "unlikely")), return(c("downstream", "unlikely")))
    }else{
      ifelse(strand(junctions_gr[i]) == "+", return(c("downstream", "unlikely")), return(c("upstream", "unlikely")))
    }
    #}
  }
  #print("inconclusive")
  return(c("inconclusive", "inconclusive"))
  #}
  ## if both lists == 0 return intergenic else inconclusive
  #return(c("inconclusive", "inconclusive"))
}

