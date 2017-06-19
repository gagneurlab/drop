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
    # Choose best among gene duplicates
    #
    proteome_data_table <- filter_proteome_for_duplicated_genes(
        proteome_data_table
    )
    print(dim(proteome_data_table))
    
    
    
    
########################################
    column_ids=c(column_protein_id, column_gene_id)
    
    # select group representative
    for(i in column_ids){
        proteome_data_table[,eval(paste0(i, "_FIRST")):=tstrsplit(
                get(i), ';', fixed=T, keep=1
            )
        ]
    }
    head(proteome_data_table)
    
    
    
    
    

    
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
    # Treat IDs
    #
    message('Number of entries in table: ', nrow(proteinGroups_data_table))
    
    noid = proteinGroups_data_table[ is.na(Protein.IDs)]
    if(nrow(noid)>0)
        message('Number of proteins without protein ID: ', nrow(noid))
    
    #
    # first_protein as KEY
    setkey(proteinGroups_data_table, first_protein)
    
    # replace 0 in intensity columns
    cols_intensity = grep(intensity_column_pattern, colnames(proteinGroups_data_table), value=T)
    stopifnot(length(cols_intensity)>0)
    replace_value_in_dt_columns(proteinGroups_data_table, 0, cols_intensity, NA)
    
    
    
    #
    # OUTLIER SAMPLES
    #
    
    # get NA frequency per sample
    sample_na_freq = sapply(res_dt[,cols_intensity, with=F], function(j) 
            sum(is.na(j))/nrow(res_dt) 
    )
    bad_samples = names(which(sample_na_freq > max_na_frequency))
    
    if(length(bad_samples)>0)
        message('Samples with NA frequency> ',max_na_frequency,':\n', paste(bad_samples, collapse='\n'))
    
    # remove poorly measured samples
    cols_intensity = setdiff(cols_intensity, bad_samples)
    res_dt = res_dt[, setdiff(colnames(res_dt), bad_samples), with=F]
    
    #
    # OUTLIER GENES
    #
    
    # remove genes with no expression
    genes_summed_expr = apply(res_dt[,cols_intensity, with=F], 1, function(i) sum(i, na.rm=T) )
    bad_genes = res_dt[which(genes_summed_expr==0), first_gene]
    if(length(bad_genes)>0)
        message('Number of proteins without expression ', length(bad_genes))
    
    res_dt = res_dt[!first_gene %in% bad_genes]
    
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







