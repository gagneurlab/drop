### 20210427 karoline lutz
### calculate splice type for aberrant splicing junctions and frameshift
### each junction in the aberrant splicing results table is compared to the most freq used intron at the junction position
### most frequently used intron is based on median count
### based on the intron the junction is compared to, it will be labelled as exonTruncation, exonSkipping, upstream/downstream,...
### theta juctions are labelled as exonic, intronic, spliceSite or up/downstream



# the main function to calculate the aberrant splice types for all junctions in the drop results table
# exonSkipping, exonTruncation, intergenic,...
aberrantSpliceType<- function(junctions_dt, fds, txdb){
  print("preparing..")
  # the first part is only done for the psi junctions, theta junctions are annotated differently
  psi_positions <- which(junctions_dt$type != "theta")
  # to create a granges object out of the dt, the STRAND column has to be renamed or else the function will get confused about having two 'strand' columns
  colnames(junctions_dt)[which(names(junctions_dt) == "STRAND")] <- "strand2"
  junctions_gr <- makeGRangesFromDataFrame(junctions_dt[psi_positions], keep.extra.columns = T)
  seqlevelsStyle(txdb) <- seqlevelsStyle(junctions_gr)
  seqlevelsStyle(fds) <- seqlevelsStyle(junctions_gr)
  exons <- exons(txdb)
  
  # get all introns that are actually expressed and contained in the fds
  # add a median count for each intron
  introns_tmp <- unique(unlist(intronsByTranscript(txdb)))
  fds_known <- fds[unique(to(findOverlaps(introns_tmp, rowRanges(fds, type = "j"), type = "equal"))),]
  grIntrons <- rowRanges(fds_known, type="j")
  introns <- as.data.table(grIntrons)
  introns <- introns[,.(seqnames, start, end, strand)]
  
  sampleCounts <- K(fds_known, type = "j")
  introns[, "medianCount" := rowMedians(as.matrix(sampleCounts))]
  intron_ranges <- makeGRangesFromDataFrame(introns, keep.extra.columns = TRUE)
  
  ## prepare the results column
  junctions_dt[, aberrantSpliceType := "NA"]
  junctions_dt[, causesFrameshift := "NA"]
  # junctions that have the same start and end as the most frequently used intron at that position
  junctions_dt[which(junctions_dt$annotatedJunction=="both"), aberrantSpliceType := "sameIntron"]
  junctions_dt[which(junctions_dt$annotatedJunction=="both"), causesFrameshift := "unlikely"]
  
  starts <- which(junctions_dt[psi_positions]$annotatedJunction=="start")
  ends <- which(junctions_dt[psi_positions]$annotatedJunction=="end")
  nones <- which(junctions_dt[psi_positions]$annotatedJunction=="none")
  
  
  print("calculating aberrant splice types")
  print("start junctions")
  
  # calculate the intron overlap first
  start_junctions <- findOverlaps(junctions_gr, intron_ranges, type = "start")
  start_results <- sapply(starts, function(i){
    # get indices of all overlaps for the current junction
    overlap <- to(start_junctions)[which(from(start_junctions) == i)]
    # get the intron with the highest median 
    maxExpr <- overlap[which.max(intron_ranges[overlap, ]$medianCount)]
    
    # start the actual comparison
    return(compareEnds(junctions_gr, i, maxExpr, F, intron_ranges, exons))
  })
  
  junctions_dt[psi_positions[starts], causesFrameshift:=start_results[2,]]
  junctions_dt[psi_positions[starts], aberrantSpliceType := start_results[1,]]
  
  
  print("end junctions")
  # calculate overlap for all junctions in this category
  end_junctions <- findOverlaps(junctions_gr, intron_ranges, type = "end")
  end_results <- sapply(ends, function(i){
    # get indices of all overlaps for the current junction
    overlap <- to(end_junctions)[which(from(end_junctions) == i)]
    # get intron with highest median
    maxExpr <- overlap[which.max(intron_ranges[overlap, ]$medianCount)]
    # start comparison
    return(compareStarts(junctions_gr, i, maxExpr, F, intron_ranges, exons))
  })
  
  junctions_dt[psi_positions[ends], causesFrameshift:=end_results[2,]]
  junctions_dt[psi_positions[ends], aberrantSpliceType := end_results[1,]]
  
  
  print("none junctions pt1")
  # calculate overlap with introns for junctions in this category
  none_junctions <- findOverlaps(junctions_gr, intron_ranges)
  none_results <- sapply(nones, function(i){
    # get indices of all overlapping introns
    overlap <- to(none_junctions)[which(from(none_junctions) == i)]
    # if overlap is 0, theres no overlaping intron -> continue with those later
    if(length(overlap) == 0) return(c("noOverlap", "inconclusive"))
    # get max expressed intron
    maxExpr <- overlap[which.max(intron_ranges[overlap, ]$medianCount)]
    # calculate type on both ends of junction
    st = compareStarts(junctions_gr, i, maxExpr, T, intron_ranges, exons)
    en = compareEnds(junctions_gr, i, maxExpr, T, intron_ranges, exons)
    
    # merge results for both ends
    # if one end is a single Exon Skipping and the other end an exon elongation -> final result is single Exon Skipping
    if((st[1] == "exonSkipping" & en[1] == "intronTruncation") ||
       (en[1] == "exonSkipping" & st[1] == "intronTruncation")){
      ## only one is single exonSkipping, the other is trunc/elong
      if((as.integer(st[3])+as.integer(en[3]))%%3 != 0){
        frs = "likely"
      }else{ frs = "unlikely"}
      return(c("exonSkipping", frs))
    } 
    
    # at least one end is labelled as exon skipping
    if(st[1] %in% c("multipleExonSkipping", "exonSkipping") || en[1] %in% c("multipleExonSkipping", "exonSkipping")) return(c("multipleExonSkipping", "inconclusive"))
    
    # no exon skipping but just normal truncation or elongation
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
  # continue with junctions without intron overlap
  # calculate their exon overlaps
  exon_junctions <- findOverlaps(junctions_gr, exons)
  noLaps_results <- sapply(noLaps, function(i){
    overlap <- to(exon_junctions)[which(from(exon_junctions) == i)]
    # if also no exon overlap -> intergenic or inconclusive
    if(length(overlap) == 0) return(checkIntergenic(junctions_gr, i, refseq.genes)) 
    
    # for the exons, check if splice site is contained in the exon
    # exon truncation only if splice site is fully contained in exon
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
    
    return(c("inconclusive","inconclusive"))
  })
  junctions_dt[psi_positions[noLaps], causesFrameshift:=noLaps_results[2,]]
  junctions_dt[psi_positions[noLaps], aberrantSpliceType := noLaps_results[1,]]
  
  
  
  # theta annotations
  print("Annotate theta junctions")
  thetas <- which(junctions_dt$type == "theta")
  junctions_gr <- makeGRangesFromDataFrame(junctions_dt[thetas,], keep.extra.columns = TRUE)
  
  # label all as intronic first if they have any intron overlap
  intronic <- unique(from(findOverlaps(junctions_gr, introns_tmp)))
  junctions_dt[thetas, aberrantSpliceType := "inconclusive"]
  junctions_dt[thetas[intronic], aberrantSpliceType := "intronic"]
  
  # for exonic, check if theta is fully contained in an exon
  # if one end is in an intron and the other in an exon it is a splice site
  print("exonic")
  exonic <- unique(from(findOverlaps(junctions_gr, exons)))
  within <- findOverlaps(junctions_gr, exons, type = "within")
  all <- findOverlaps(junctions_gr, exons)
  exonic_results <- sapply(exonic, function(i){
    w <- unique(to(within)[which(from(within) == i)])
    a <- unique(to(all)[which(from(all) == i)])
    if(length(a) == length(w)) return("exonic")
    return("spliceSite")
  })
  
  junctions_dt[thetas[exonic], aberrantSpliceType := exonic_results]
  
  # check cases that don't overlap with an exon/intron
  print("nones")
  nones <- which(junctions_dt[thetas,]$aberrantSpliceType == "inconclusive")
  none_results <- sapply(nones, function(i){
    if(length(findOverlaps(junctions_gr[i], refseq.genes)) > 0) return("inconclusive")
    #return("intergenic")
    if(start(refseq.genes[nearest(junctions_gr[i], refseq.genes)]) > start(junctions_gr[i])){
      ifelse(strand(junctions_gr[i]) == "+", return("upstream"), return("downstream"))
    }else{
      ifelse(strand(junctions_gr[i]) == "+", return("downstream"), return("upstream"))
    }
  })
  junctions_dt[thetas[nones], aberrantSpliceType := none_results]
  print("thetas done")
  
  # add distance to closest neighbour gene for intergenic results (both psi and theta)
  print("adding distances to nearest gene")
  up <- which(junctions_dt$aberrantSpliceType == "upstream")
  down <- which(junctions_dt$aberrantSpliceType == "downstream")
  
  # create full grange object containing psi and theta
  junctions_gr <- makeGRangesFromDataFrame(junctions_dt, keep.extra.columns = T)
  
  print("Calculate distances")
  junctions_dt[, distNearestGene := 0]
  if(length(up) > 0){
    distanceNearestGene_up <- sapply(up, function(i) min(distance(junctions_gr[i], refseq.genes), na.rm = T))
    if(length(distanceNearestGene_up > 0)){
      junctions_dt[up, distNearestGene := distanceNearestGene_up]
    } else{
      junctions_dt[up, distNearestGene := NA]
      print("No distances found for upstream")
    }
  }else{print("No upstream targets")}
  
  if(length(down) > 0){
    distanceNearestGene_down <- sapply(down, function(i) min(distance(junctions_gr[i], refseq.genes), na.rm = T))
    if(length(distanceNearestGene_down > 0)){
      junctions_dt[down, distNearestGene := distanceNearestGene_down]
    }else{
      junctions_dt[down, distNearestGene := NA]
      print("No distances found for downstream")
    }
  }else{print("No downstream targets")}
  
  # change column name back to original
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
    
    if(((ss_start - intron_start)%%3) != 0){
      frs = "likely"
    }else{ frs = "unlikely"}
    
    ifelse(shift_needed, return(c("intronTruncation", frs, (ss_start - intron_start))), return(c("intronTruncation", frs)))
    
  }
  
  ## check if splice site ends in following exon -> exon truncation -> FRS -> done
  if(intron_start > ss_start){
    
    ## create dummy exon find all exons starting from that intron end
    dummy_exon <- GRanges(
      seqnames = toString(seqnames(intron_ranges[max_lap])),
      ranges = IRanges(intron_start-2, end = intron_start -1),
      strand = toString(strand(intron_ranges[max_lap]))
    )
    exonChoices <- to(findOverlaps(dummy_exon, exons, type = "end"))
    
    for(j in exonChoices){
      exon_start = start(exons[j])
      
      if(exon_start <= ss_start){
        if((end(exons[j]) - ss_start + 1)%%3 != 0){
          frs = "likely"
        }else{frs = "unlikely"}
        ifelse(shift_needed, return(c("exonTruncation", frs, (-1)*(end(exons[j]) - ss_start + 1))), return(c("exonTruncation", frs)))
        
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
        ifelse(shift_needed, return(c("multipleExonSkipping", "inconclusive", 0)), return(c("multipleExonSkipping", "inconclusive")))
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
        ifelse(shift_needed, return(c("multipleExonSkipping", "inconclusive", 0)), return(c("multipleExonSkipping", "inconclusive")))
      }
      
      if(ss_start >= start(intron_ranges[secItrChoices[maxExpr]])){
        ## check if there is no other exon in that range
        if(length(findOverlaps(exons, intron_ranges[secItrChoices[maxExpr]], type = "within")) == 0){
          ## clear exon skipping, only exon is skipped 
          ## calculate frameshift, skipped exon plus possible exon elongation 
          
          shift = (-1)*(end(exons[exonChoices[1]]) - start(exons[exonChoices[1]]) + 1) + ss_start - start(intron_ranges[secItrChoices[maxExpr]])
          
          frs = ifelse(shift%%3 == 0,"unlikely","likely")
          ifelse(shift_needed, return(c("exonSkipping", "inconclusive", shift)), return(c("exonSkipping", frs)))
        }
      }
    } ## single exon skipping end
    
  }
  
  ## splice site longer than one intron + exon -> not defined for now
  ifelse(shift_needed, return(c("multipleExonSkipping", "inconclusive", 0)), return(c("multipleExonSkipping", "inconclusive")))
  
}


compareEnds <- function(junctions_gr, i, max_lap, shift_needed, intron_ranges, exons){
  intron_end = end(intron_ranges[max_lap])
  ss_end = end(junctions_gr[i])
  
  ## found the most freq intron with same start again
  ## check if intron ends after splice site -> exon elongation -> FRS -> done
  if(intron_end > ss_end){
    
    if(((intron_end - ss_end)%%3) != 0){
      frs = "likely"
    }else{ frs = "unlikely"}
    
    ifelse(shift_needed, return(c("intronTruncation", frs, (intron_end - ss_end))), return(c("intronTruncation", frs)))
    
  }
  
  ## check if splice site ends in following exon -> exon truncation -> FRS -> done
  if(intron_end < ss_end){
    
    ## create dummy exon find all exons starting from that intron end
    dummy_exon <- GRanges(
      seqnames = toString(intron_ranges[max_lap]@seqnames@values),
      ranges = IRanges(intron_end + 1, end = intron_end + 2),
      strand = toString(intron_ranges[max_lap]@strand@values)
    )
    exonChoices <- to(findOverlaps(dummy_exon, exons, type = "start"))
    
    for(j in exonChoices){
      exon_end = end(exons[j])
      
      if(exon_end >= ss_end){
        if((ss_end - start(exons[j]) + 1)%%3 != 0){
          frs = "likely"
        }else{frs = "unlikely"}
        ifelse(shift_needed, return(c("exonTruncation",frs, (-1)*(ss_end - start(exons[j]) + 1))), return(c("exonTruncation",frs)))
        
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
        ifelse(shift_needed, return(c("multipleExonSkipping", "inconclusive", 0)), return(c("multipleExonSkipping", "inconclusive")))
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
        ifelse(shift_needed, return(c("multipleExonSkipping", "inconclusive", 0)), return(c("multipleExonSkipping", "inconclusive")))
      }
      
      if(ss_end <= end(intron_ranges[secItrChoices[maxExpr]])){
        ## check if there is no other exon in that range
        if(length(findOverlaps(exons, intron_ranges[secItrChoices[maxExpr]], type = "within")) == 0){
          ## clear exon skipping, only exon is skipped 
          ## calculate frameshift, skipped exon plus possible exon elongation at end
          shift = (-1)*(end(exons[exonChoices[1]]) - start(exons[exonChoices[1]]) + 1) + end(intron_ranges[secItrChoices[maxExpr]]) - ss_end
          frs = ifelse(shift%%3 == 0,"unlikely","likely")
          ifelse(shift_needed, return(c("exonSkipping", "inconclusive", shift)), return(c("exonSkipping", frs)))
        }
      }
    } ## single exon skipping end
    
    
  }
  
  ## splice site longer than one intron + exon -> not defined for now
  ifelse(shift_needed, return(c("multipleExonSkipping", "inconclusive", 0)), return(c("multipleExonSkipping", "inconclusive")))
  
}


checkIntergenic <- function(junctions_gr, i, refseq.genes){
  ## check if start > 1000
  ## start - 1000, end + 1000
  start = start(junctions_gr[i]) 
  
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
  
  dist = min(distance(test_junction, refseq.genes), na.rm = T)
  
  if(dist > 0){
    if(start(refseq.genes[nearest(junctions_gr[i], refseq.genes)]) > start){
      ifelse(strand(junctions_gr[i]) == "+", return(c("upstream", "unlikely")), return(c("downstream", "unlikely")))
    }else{
      ifelse(strand(junctions_gr[i]) == "+", return(c("downstream", "unlikely")), return(c("upstream", "unlikely")))
    }
    
  }
  
  return(c("inconclusive", "inconclusive"))
  
}

