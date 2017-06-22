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
read_data_table_from_proteinGroupsTxt <- function(
    inputfile,
    intensity_column_pattern='LFQ.intensity.',
    column_ids=c('Protein.IDs', "Gene.names"),
    column_sample_id='PROTEOME_ID',
    column_intensity='LFQ_INTENSITY',
    replace_zero_with_NA=TRUE
){
    stopifnot(exists('replace_value_in_dt_columns'))
    library(tidyr)
    library(data.table)
    
    #
    # read file
    #
    
    # use read.delim for big integers
    proteinGroups_data_table= data.table(read.delim( inputfile, 
            comment.char='#', na.strings=c('',NaN, 'NA', 'NaN'), stringsAsFactors = FALSE, as.is=T
        ))
    
    # remove multi dots
    setnames( proteinGroups_data_table, gsub('\\.\\.+', '\\.', colnames(proteinGroups_data_table)))
    # remove trailing dots
    setnames( proteinGroups_data_table, gsub('\\.+$', '', colnames(proteinGroups_data_table)))
    
    
    #
    # process table
    #
    
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
    if(length(cols_intensity)==0){
        stop(paste("Bad intensity search pattern. No columns selected.\n",
                colnames(proteinGroups_data_table))
        )
    }
    
    
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
    pdt_tidy <- as.data.table(gather_(
            proteinGroups_data_table, 
            key=column_sample_id, 
            value=column_intensity, 
            gather_cols = cols_intensity
        )
    )
    
    # remove intensity prefix
    pdt_tidy[
        ,c(column_sample_id):=list(
            gsub(intensity_column_pattern, '', get(column_sample_id))
        )]
    
    # rename ID columns: remove plural s, replace '.' by '_', make upper case
    for(i in column_ids){
        setnames(pdt_tidy, 
            sub(i, 
                toupper(sub('s$', '', sub('\\.', '_', i))), 
                names(pdt_tidy)
            )
        )
    }
    
    # return
    return(pdt_tidy)
}







