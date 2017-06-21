# R function
# 
# Author: baderd
###############################################################################


get_unique_protein_ids <- function(
    proteome_data_table,
    column_protein_id='PROTEIN_ID'
){
    unique(proteome_data_table[, get(column_protein_id)])
}


get_unique_gene_ids <- function(
    proteome_data_table,
    column_genen_id='GENE_NAME'
){
    unique(proteome_data_table[, get(column_gene_id)])
}



#' proteins with no gene name 
#' 
#' @return data.table; subset of input where column_gene_id is NA
#' 
get_proteome_with_unknown_genes <- function(
    proteome_data_table,
    column_protein_id='PROTEIN_ID',
    column_gene_id='GENE_NAME'
){
    res <- proteome_data_table[
        is.na(get(column_gene_id)), 
        c(column_protein_id, column_gene_id), 
        with=F
    ]
    return(unique(res))
}



#' proteins with duplicated gene names
#' 
get_proteome_with_duplicated_genes <- function(
    proteome_data_table,
    column_protein_id='PROTEIN_ID',
    column_gene_id='GENE_NAME'
){
    # get unique matching protein to gene
    dt <- unique(
        proteome_data_table[,c(column_protein_id, column_gene_id), with=F]
    )
    index_duplicated <- duplicated(dt[, get(column_gene_id)])
    
    # names of duplicates
    genes_duplicated <- get_unique_gene_ids(dt[index_duplicated, ])
    
    # return all affected entries
    res <- proteome_data_table[get(column_gene_id) %in% genes_duplicated]
    return(res)
}



#' select_best_among_duplicates
#' 
#' For all protein entries of one gene, pick best. 
#' 
select_best_among_duplicates <- function(
    duplicated_proteome_dt,
    column_protein_id='PROTEIN_ID',
    column_gene_id='GENE_NAME',
    column_intensity='LFQ_INTENSITY'
){
    # ensure that only 1 gene is compared
    num_genes_tested <- length(unique(duplicated_proteome_dt[,get(column_gene_id)]))
    stopifnot(num_genes_tested==1)
    
    # mean intensity by protein across samples
    info_by_pid <- duplicated_proteome_dt[, 
        mean(get(column_intensity), na.rm=TRUE), 
        by=eval(column_protein_id)
    ]
    
    # return highest mean
    best_pid <- info_by_pid[V1==max(V1), get(column_protein_id)]
    res <- duplicated_proteome_dt[get(column_protein_id)==best_pid]
    res
}


#' filter_proteome_for_duplicated_genes
#' 
#' Return subset of input proteome. 
#' Only best genes among duplicates are returned. 
#' 
filter_proteome_for_duplicated_genes <- function(
    proteome_data_table,
    column_protein_id='PROTEIN_ID',
    column_gene_id='GENE_NAME',
    column_intensity='LFQ_INTENSITY'
){
    dupl_gene_pid <- get_unique_protein_ids(
        get_proteome_with_duplicated_genes(
            proteome_data_table, column_protein_id, column_gene_id
        )
    )
    dupl_gene_names <- get_unique_gene_ids(proteome_data_table[
            get(column_protein_id) %in% dupl_gene_pid
        ]
    )
    # report duplicates lost
    num_duplicates_lost = length(dupl_gene_pid) - length(dupl_gene_names)
    message('Number of duplicates gene names: ', length(dupl_gene_names))
    message('Number of protein duplicates removed: ', num_duplicates_lost)
    
    
    # select best for each duplicated gene
    best_dupl_proteome_dt <- data.table()
    for(g in dupl_gene_names){
        dupl_dt <- proteome_data_table[get(column_gene_id)==g]
        res <- select_best_among_duplicates(
            dupl_dt, column_protein_id, column_gene_id, column_intensity
        )
        best_dupl_proteome_dt <- rbind(best_dupl_proteome_dt, res)
    }
    
    # remove all duplicate entries
    proteome_data_table <- proteome_data_table[
        !get(column_protein_id) %in% dupl_gene_pid
    ]
    # add best duplicates back 
    proteome_data_table <- rbind(proteome_data_table, best_dupl_proteome_dt)
}


replace_groupid_with_firstid <- function(
    proteome_data_table,
    columns_groupid=c('PROTEIN_ID', 'GENE_NAME'),
    remove_groupid_columns=TRUE,
    separator=';'
){
    library(data.table)
    first_suffix <- "_FIRST"
    group_prefix <- 'GROUP_'
    
    # select group representative
    for(i in columns_groupid){
        proteome_data_table[,eval(paste0(i, first_suffix)):=tstrsplit(
                get(i), separator, fixed=T, keep=1
            )
        ]
    }
    
    # rename columns
    for(i in columns_groupid){
        setnames(
            proteome_data_table, 
            sub(paste0(i,'$'), paste0(group_prefix, i), names(proteome_data_table))
        )
    }
    setnames(
        proteome_data_table, 
        gsub(first_suffix, "", names(proteome_data_table))
    )
    
    # return with our without group ids
    if(remove_groupid_columns){
        # remove GROUP_
        proteome_data_table[,grep(group_prefix, names(proteome_data_table), value=T):=list(NULL)]
    }
    return(proteome_data_table)
}

#' filter_proteome_for_genes
#' 
filter_proteome_by_gene_properties <- function(
    proteome_data_table
){
    column_protein_id='PROTEIN_ID'
    column_gene_id='GENE_NAME'
    column_intensity='LFQ_INTENSITY'
    remove_proteins_without_gene_name=TRUE
    
    
    print(dim(proteome_data_table))
    
    # 
    # Remove proteins w/o gene name
    #
    if(remove_proteins_without_gene_name){
        nogene_pid <- get_unique_protein_ids(get_proteome_with_unknown_genes(
            proteome_data_table, column_protein_id, column_gene_id
        ))
        proteome_data_table <- proteome_data_table[
            !get(column_protein_id) %in% nogene_pid
        ]
        message("Number of proteins without gene name: ", length(nogene_pid))
    }
    print(dim(proteome_data_table))
    
    
    #
    # Keep only first ID of protein group
    #
    proteome_data_table <- replace_groupid_with_firstid(proteome_data_table)
    
    
    #
    # Choose best among gene duplicates
    #
    proteome_data_table <- filter_proteome_for_duplicated_genes(
        proteome_data_table
    )
    print(dim(proteome_data_table))
}



compute_sample_na_frequency <- function(
    proteome_data_table,
    column_intensity='LFQ_INTENSITY',
    column_sample='PROTEOME_ID'
){
    proteome_data_table[,
        SAMPLE_NA_FREQ := sum(is.na(get(column_intensity)))/.N, 
        by=get(column_sample)
    ]
}


compute_gene_na_frequency <- function(
    proteome_data_table,
    column_intensity='LFQ_INTENSITY',
    column_id='GENE_NAME'
){
    proteome_data_table[,
        GENE_NA_FREQ := sum(is.na(get(column_intensity)))/.N, 
        by=get(column_id)
    ]
}



#' proteinGroups_quality_control
#' 
#' @param proteinGroups_data_table result of proteinGroups_read_table
#' 
proteinGroups_dt_quality_control = function( proteinGroups_data_table, 
    intensity_column_pattern= 'LFQ', 
    max_na_frequency= 0.5, 
    low_expr_quantile=NULL, 
    return_matrix=TRUE, 
    verbose=FALSE
){
   
    #
    # OUTLIER GENES
    #
    

    # remove genes with low expression (treat NA as 0)
    if(!is.null(low_expr_quantile)){
        ugly_idx = get_low_expression_indices(
            replace_na(res_dt[,cols_intensity, with=F],0), low_expr_quantile, 1
        )
        if(sum(ugly_idx)>0)
            message('Number of proteins with ',low_expr_quantile,' quantile not measured: ', sum(ugly_idx))
        
        res_dt= res_dt[!ugly_idx]
    }else{
        message('low_expr_quantile=',low_expr_quantile,'. Do you want to keep proteins with very low expression?')
    }
    
    #
    # RETURN filtered result
    #
    message('Final number of quality approved genes/proteins: ', nrow(res_dt))
    
    # convert to gene x sample matrix?
    if(return_matrix){
        m= as.matrix(res_dt[,cols_intensity, with=F])
        rownames(m)= res_dt[,first_gene]
        return(m)
    }else
        return(res_dt)
}







