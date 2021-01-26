
# Function to add HPO terms to results table

add_HPO_cols <- function(RES, sample_id_col = 'sampleID', 
                         gene_name_col = 'hgncSymbol', hpo_file = NULL){
  require(data.table)
  
  filename <- ifelse(is.null(hpo_file), 
                 'https://www.cmm.in.tum.de/public/paper/drop_analysis/resource/hpo_genes.tsv.gz',
                 hpo_file)
  
  hpo_dt <- fread(filename)
  
  # change column names
  setnames(RES, old = c(sample_id_col, gene_name_col), new = c('sampleID', 'hgncSymbol'))
  
  f2 <- merge(RES[, .(sampleID, hgncSymbol)], 
            hpo_dt[,.(hgncSymbol, HPO_id, HPO_label)], by = 'hgncSymbol')
  
  if(nrow(f2) > 0){
    f3 <- merge(f2, sa[,.(RNA_ID, HPO_TERMS)], by.x = 'sampleID', by.y = 'RNA_ID')
    f3[, HPO_match := HPO_id %in% unlist(strsplit(HPO_TERMS, split = ','))]
    f3 <- f3[HPO_match == TRUE]
    if(nrow(f3) > 0){
      f4 <- f3[, .(HPO_id_overlap = paste(HPO_id, collapse = ', '), 
                 HPO_label_overlap = paste(HPO_label, collapse = ', ')), 
             by = .(sampleID, hgncSymbol)]
      RES <- merge(RES, f4, by = c('sampleID', 'hgncSymbol'), all.x = TRUE)
    }
  }
  
  setnames(RES, old = c('sampleID', 'hgncSymbol'), new = c(sample_id_col, gene_name_col))
  return(RES)
}
