# R function
#
# Author: baderd
###############################################################################


summarize_columns_of_matrix  <- function(
    ematrix, old2new_column_matching
){
    # remove old samples without match in new columns
    new_na = which(is.na(old2new_column_matching))
    if(length(new_na)>0)
        message('No match for these columns: ', names(new_na))
    old2new_column_matching = na.omit(old2new_column_matching)
    
    # count replicates to be summarized
    tmp_new_counts_per_old_columns= table(old2new_column_matching)
    
    ematrix_new_columns = sapply(unique(old2new_column_matching), function(old_coln){
            if(tmp_new_counts_per_old_columns[old_coln]==1){
                ematrix[, names(which(old2new_column_matching==old_coln))]
            }else 
            # compute mean over multiple columns in log space
            if(tmp_new_counts_per_old_columns[old_coln] >1){
                2^rowMeans( log2(ematrix[, names(which(old2new_column_matching==old_coln))]), na.rm=T)
            }else{
                stop('Bad new name')
            }
        })
    rownames(ematrix_new_columns)= rownames(ematrix)
    return(ematrix_new_columns)
}


#' summarize_ematrix_by_fibroblast
#' 
#' For specific application to convert rna/protein matrices to fibroblast matrices
#' 
summarize_ematrix_by_fibroblast = function( ematrix, converter_fct_name, ...){
	summarize_columns_of_matrix(
	    ematrix, 
	    sapply(colnames(ematrix), get(converter_fct_name),...)
    )
}
