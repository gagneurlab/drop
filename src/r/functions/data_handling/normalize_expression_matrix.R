

#' normalize_expression_matrix
#' 
#' @param ematrix numeric matrix, containing expression values
#' @param sizefactor logical (default TRUE); should DESeq2 size factor normalization be applied?
#' @param rowcenter logical (default TRUE); should the matrix be normalized with respect to mean per row?
#' @param log2scale logical (default TRUE); should the matrix be transformed to log2 space? 
#' 		if it is already in log space, set FALSE, sizefactor and rowcenter will still work
#' @param nonzero numeric (default 1); to avoid log2(0)=-Inf, especially for RNA counts
#' 
normalize_expression_matrix= function(ematrix, sizefactor=TRUE, rowcenter=TRUE, log2scale=TRUE, nonzero=1
){
	require(DESeq2)
	
	if(sizefactor){
		sf= DESeq2::estimateSizeFactorsForMatrix(ematrix)
		ematrix= t(t(ematrix)/sf)
	}
	
	if(log2scale){
		ematrix= log2(nonzero + ematrix)
	}
	
	if(rowcenter){
		if(log2scale){
			ematrix= ematrix - rowMeans(ematrix, na.rm=T)
		}else{
			if(any(ematrix <0, na.rm=T)){
				ematrix= ematrix - rowMeans(ematrix, na.rm=T)
			}else{
				ematrix= ematrix / rowMeans(ematrix, na.rm=T)
			}
		}
	}
	return(ematrix)
}


normalize_ematrix_by_feature <- function(ematrix, anno_table, feature_name, key_anno
){
    require(data.table)
    setkeyv(anno_table, key_anno)
    norm_vector = anno_table[rownames(ematrix), get(feature_name)]
    norm_mat = ematrix / norm_vector
    return(norm_mat)
}
# normalize_ematrix_by_feature(rna_fibro_lvl, ANNO_ALMANAC, 'exon_length', 'Gene.names')



