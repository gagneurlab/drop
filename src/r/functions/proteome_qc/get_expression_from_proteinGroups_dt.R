#'---
#' title: R script
#' author: Daniel Bader
#' output: 
#'   html_document:
#'     toc_float: yes
#'     toc: yes
#'---
#'



#' Get expression table from proteinGroups.txt based data.table.
#' 
#' @param column_ids character vector; first entry should point to unique IDs
#' @return data.table
#' 
get_expression_from_proteinGroups_dt <- function(
    proteinGroups_data_table,
    intensity_column_pattern='LFQ.intensity.',
    column_ids=c('Protein.IDs', "Gene.names"),
    column_sample_id='PROTEOME_ID',
    column_intensity='LFQ_INTENSITY',
    replace_zero_with_NA=TRUE
){
    stopifnot(exists('replace_value_in_dt_columns'))
    library(tidyr)
    library(data.table)
    
    # remove rows without ID
    noid <- is.na(proteinGroups_data_table[, get(column_ids[1])])
    if(any(noid)){
        message('Number of proteins without protein ID: ', nrow(noid))
    }
    proteinGroups_data_table <- proteinGroups_data_table[!noid]
    
    # grep intensity columns
    cols_intensity = grep(
        intensity_column_pattern, 
        colnames(proteinGroups_data_table), 
        value=TRUE
    )
    stopifnot(length(cols_intensity)>0)
    
    # reduce columns
    proteinGroups_data_table <- proteinGroups_data_table[ 
        ,c(column_ids, cols_intensity), with=F
    ]
    
    # insert NA? replace 0 in intensity columns
    if(replace_zero_with_NA){
        replace_value_in_dt_columns(
            proteinGroups_data_table, 0, cols_intensity, NA
        )
    }
    
    # make tidy expression table
    protein_intensity_dt <- as.data.table(gather_(
            proteinGroups_data_table, 
            key=column_sample_id, 
            value=column_intensity, 
            gather_cols = cols_intensity
        )
    )
    
    # remove intensity prefix
    protein_intensity_dt[
        ,c(column_sample_id):=list(
            gsub(intensity_column_pattern, '', get(column_sample_id))
        )]
    
    # rename ID columns: remove plural s, replace '.' by '_', make upper case
    for(i in column_ids){
        setnames(protein_intensity_dt, 
            sub(i, 
                toupper(sub('s$', '', sub('\\.', '_', i))), 
                names(protein_intensity_dt)
            )
        )
    }
    
    # return
    protein_intensity_dt
}







