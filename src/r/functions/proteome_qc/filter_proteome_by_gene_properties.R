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
    column_gene_id='GENE_NAME'
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





#' filter_proteome_for_genes
#' 
#' meta function to call multiple sub filters
#' 
filter_proteome_by_gene_properties <- function(
    proteome_data_table,
    column_protein_id='PROTEIN_ID',
    column_gene_id='GENE_NAME',
    column_intensity='LFQ_INTENSITY',
    remove_proteins_without_gene_name=TRUE,
    keep_first_id_only=TRUE,
    verbose=FALSE
){
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
        if(verbose){
            message("Number of proteins without gene name: ", length(nogene_pid))
        }
    }
    
    #
    # Keep only first ID of protein group
    #
    if(keep_first_id_only){
        proteome_data_table <- replace_groupid_with_firstid(
            proteome_data_table,
            columns_groupid = c(column_protein_id, column_gene_id),
            remove_groupid_columns = TRUE,
            separator = ';'
        )
    }
    
    
    #
    # Choose best among gene duplicates
    #
    proteome_data_table <- filter_proteome_for_duplicated_genes(
        proteome_data_table,
        column_protein_id = column_protein_id,
        column_gene_id = column_gene_id,
        column_intensity = column_intensity
    )
    
    # END
    return(proteome_data_table)
}








