
#' 
#' # Convert to matrix
#' 
convert_tidy_table_to_numeric_matrix <- function(
    tidytable, 
    coln_rownames, 
    coln_colnames, 
    coln_values
){
    library(data.table)
    mat <- dcast(
        data = tidytable, 
        formula = as.formula(paste(coln_rownames, '~', coln_colnames)),
        value.var = coln_values
    )
    rn <- mat[,get(coln_rownames)]
    mat[,eval(coln_rownames):=NULL]
    mat <- as.matrix(mat)
    rownames(mat) <- rn
    return(mat)
}

# pmat <- convert_tidy_table_to_numeric_matrix(pdt, 'GENE_NAME', 'PROTEOME_ID', 'LFQ_INTENSITY')
# head(pmat)
