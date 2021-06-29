### 20210525 karoline lutz
### validate the exonSkippings non overlap cases, check for beyond gene/multigenic/intergenic

checkExonSkipping <- function(junctions_dt, txdb){
  psi_positions <- which(junctions_dt$type != "theta")
  colnames(junctions_dt)[which(names(junctions_dt) == "STRAND")] <- "strand2"
  junctions_gr <- makeGRangesFromDataFrame(junctions_dt[psi_positions], keep.extra.columns = T)
  seqlevelsStyle(txdb) <- seqlevelsStyle(junctions_gr)
  
  refseq.genes<- genes(txdb)
  
  exonSkip <- which(junctions_dt[psi_positions]$aberrantSpliceType %in% c("exonSkipping", "singleExonSkipping"))
  
  print("start checking exonSkipping")
  newSkip_results <- sapply(exonSkip, function(i){
    #print(i)
    #print("----------------------------------------")
    start = start(junctions_gr[i]) 
    end = end(junctions_gr[i])
    
    #reduce the junction width so adjacent genes have a distance of 1
    if(start + 2 < end){
      start = start + 1
      end = end - 1
    }
    
    test_start <- GRanges(
      seqnames = seqnames(junctions_gr[i]),
      strand = strand(junctions_gr[i]),
      ranges = IRanges(start, start + 1)
    )
    
    test_end <- GRanges(
      seqnames = seqnames(junctions_gr[i]),
      strand = strand(junctions_gr[i]),
      ranges = IRanges(end - 1, end)
    )
    
    ## check for which genes distance to start is 0
    start_genes <- which(distance(test_start, refseq.genes) == 0)
    ## start is not in a gene
    if(length(start_genes) == 0) return("beyondGene")
    
    ## start is in a gene -> is end in same gene
    for(to in start_genes){
      ## end is in same gene
      if(distance(test_end, refseq.genes[to]) == 0) return("exonSkipping")
    }
    
    end_genes <- which(distance(test_end, refseq.genes) == 0)
    ## end is not in a gene
    if(length(end_genes) == 0) return("beyondGene")
    ## end is in a different gene
    return("multigenic")
  })
  
  print("checking exonSkipping done")
  if(length(exonSkip) > 0){
    junctions_dt[psi_positions[exonSkip], aberrantSpliceType2 := newSkip_results]
    junctions_dt[which(junctions_dt$aberrantSpliceType2 == "beyondGene"), aberrantSpliceType := "beyondGene"]
    junctions_dt[which(junctions_dt$aberrantSpliceType2 == "beyondGene"), causesFrameshift := "inconclusive"]
    junctions_dt[which(junctions_dt$aberrantSpliceType2 == "multigenic"), aberrantSpliceType := "multigenic"]
    junctions_dt[which(junctions_dt$aberrantSpliceType2 == "multigenic"), causesFrameshift := "inconclusive"]
    junctions_dt[, aberrantSpliceType2 := NULL]
  }
  
  colnames(junctions_dt)[which(names(junctions_dt) == "STRAND")] <- "strand2"
  return(junctions_dt)
}


checkInconclusive <- function(junctions_dt, txdb){
  psi_positions <- which(junctions_dt$type != "theta")
  colnames(junctions_dt)[which(names(junctions_dt) == "STRAND")] <- "strand2"
  junctions_gr <- makeGRangesFromDataFrame(junctions_dt[psi_positions], keep.extra.columns = T)
  seqlevelsStyle(txdb) <- seqlevelsStyle(junctions_gr)
  
  refseq.genes<- genes(txdb)
  
  inconclusive <- which(junctions_dt[psi_positions]$aberrantSpliceType == "inconclusive")
  print("start checking inconclusive")
  
  inconclusive_results <- sapply(inconclusive, function(i){
    #print(i)
    #print("----------------------------------------")
    start = start(junctions_gr[i]) 
    end = end(junctions_gr[i])
    
    #reduce the junction width so adjacent genes have a distance of 1
    if(start + 2 < end){
      start = start + 1
      end = end - 1
    }
    
    test_start <- GRanges(
      seqnames = seqnames(junctions_gr[i]),
      strand = strand(junctions_gr[i]),
      ranges = IRanges(start, start + 1)
    )
    
    test_end <- GRanges(
      seqnames = seqnames(junctions_gr[i]),
      strand = strand(junctions_gr[i]),
      ranges = IRanges(end - 1, end)
    )
    
    ## check for which genes distance to start is 0
    start_genes <- which(distance(test_start, refseq.genes) == 0)
    ## start is not in a gene
    if(length(start_genes) == 0) return("beyondGene")
    
    ## start is in a gene -> is end in same gene
    for(to in start_genes){
      ## end is in same gene
      if(distance(test_end, refseq.genes[to]) == 0) return("inconclusive")
    }
    
    end_genes <- which(distance(test_end, refseq.genes) == 0)
    ## end is not in a gene
    if(length(end_genes) == 0) return("beyondGene")
    ## end is in a different gene
    return("multigenic")
  })
  
  colnames(junctions_dt)[which(names(junctions_dt) == "strand2")] <- "STRAND"
  print("done checking inconclusive")
  #table(inconclusive_results)
  if(length(inconclusive) > 0) junctions_dt[psi_positions[inconclusive], aberrantSpliceType := inconclusive_results]
  
  return(junctions_dt)
}
