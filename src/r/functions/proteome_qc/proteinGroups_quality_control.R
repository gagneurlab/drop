# R function
# 
# Author: baderd
###############################################################################



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
    noid <- is.na(proteinGroups_data_table[ ,c(column_ids), with=F])
    proteinGroups_data_table <- proteinGroups_data_table[!noid]
    
    # grep intensity columns
    cols_intensity = grep(
        intensity_column_pattern, 
        colnames(proteinGroups_data_table), 
        value=TRUE
    )
    stopifnot(length(cols_intensity)>0)
    
    # insert NA? replace 0 in intensity columns
    if(replace_zero_with_NA){
        replace_value_in_dt_columns(
            proteinGroups_data_table, 0, cols_intensity, NA
        )
    }
    
    # make tidy expression table
    proteinGroups_data_table <- proteinGroups_data_table[ 
        ,c(column_ids, cols_intensity), with=F
    ]
    protein_intensity_dt <- as.data.table(gather_(
        proteinGroups_data_table, 
        key=column_sample_id, 
        value=column_intensity, 
        gather_cols = cols_intensity
    ))
    
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
    return(protein_intensity_dt)
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
    
    # set first protein and gene of group as the only ID, since they are dominating their group
    tmpid = sapply( strsplit(proteinGroups_data_table$Protein.IDs,';'), '[',1)
    proteinGroups_data_table[ ,first_protein:= tmpid ]
    tmpid = sapply( strsplit(proteinGroups_data_table$Gene.names,';'), '[',1)
    proteinGroups_data_table[ ,first_gene:= tmpid ]
    
    # first_protein as KEY
    setkey(proteinGroups_data_table, first_protein)
    
    # replace 0 in intensity columns
    cols_intensity = grep(intensity_column_pattern, colnames(proteinGroups_data_table), value=T)
    stopifnot(length(cols_intensity)>0)
    replace_value_in_dt_columns(proteinGroups_data_table, 0, cols_intensity, NA)
    
    
    #
    # Focus on genes
    #
    
    # remove proteins w/o gene name
    nogene = is.na(proteinGroups_data_table$first_gene)
    proteinGroups_data_table = proteinGroups_data_table[!nogene]
    message('Number of proteins without gene ID: ', sum(nogene))
    if(verbose)
        print(proteinGroups_data_table[nogene, first_protein])
    
    # detect duplicated gene names
    genes_duplicated = unique(proteinGroups_data_table[duplicated(first_gene), first_gene])
    
    # remove all duplicated genes
    res_dt = proteinGroups_data_table[ !first_gene %in% genes_duplicated, ]
    
    # find protein with highest avg expression and least NA values among duplicates
    best_couples = sapply( genes_duplicated, function(g){
            pids= proteinGroups_data_table[first_gene==g, first_protein]
            tmp_metrics = sapply(pids, function(p){
                    expr_single_prot = log2(as.numeric(proteinGroups_data_table[first_protein==p, cols_intensity, with=F]))
                    # compute comparables
                    c( mean= mean(expr_single_prot, na.rm=T), 
                        median= median(expr_single_prot, na.rm=T), 
                        notnafreq= 1-sum(is.na(expr_single_prot))/length(cols_intensity) 
                    )
                })
            # get protein with most maxima
            best= names(sort(colSums(tmp_metrics==rowMax(replace_na(tmp_metrics, -Inf))), decreasing=T))[1]
            return(best)
        })
    # add best proteins for duplicated genes
    res_dt = rbind(res_dt, proteinGroups_data_table[best_couples])
    setkey(res_dt, first_protein)
    stopifnot( anyDuplicated(res_dt$first_gene)==0 )
    
    # report duplicates lost
    num_duplicates_lost = nrow(proteinGroups_data_table) - nrow(res_dt)
    message('Number of protein duplicates removed: ', num_duplicates_lost, '. Keep only most abundant among duplicates.')
    
    
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







