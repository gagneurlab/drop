# R function
# 
# Author: baderd
###############################################################################


get_proteome_design_matrix = function(
    patient_fibro= NULL,
    proteome_ids = colnames(prot_intensity),
    sample_anno_columns = c('FIBROBLAST_ID', 'RNA_HOX_GROUP', 'GENDER'),
    sample_anno_dt = MITOMAP_DATATABLE,
    binarize_fibro_id= TRUE,
    return_contrasts= FALSE,
    verbose= FALSE
){
    stopifnot(!is.null(proteome_ids))
    setkey(sample_anno_dt, PROTEOME_ID)
    design_model_dt = unique(
        sample_anno_dt[ proteome_ids, c('PROTEOME_ID', sample_anno_columns), with=F ]
    )
    
    # tmp_df_rownames = design_model_dt$FIBROBLAST_ID
    if("FIBROBLAST_ID" %in% colnames(design_model_dt) && binarize_fibro_id){
        design_model_dt[,FIBROBLAST_ID:= factor(ifelse(
                FIBROBLAST_ID==patient_fibro, 'case', 'control'
            ), levels=c('control', 'case'))]
    }
    design_model_df = data.frame( design_model_dt, row.names= design_model_dt$PROTEOME_ID)
    # design_model_df = data.frame( design_model_dt, row.names= tmp_df_rownames)
    
    if(verbose){
        print(table(design_model_df[,sample_anno_columns]))
    }
    
# make formula and model matrix
    design_mat = model.matrix( 
        object= formula( paste("~",paste(sample_anno_columns, collapse="+")) ), 
        data= design_model_df
    )
    colnames(design_mat) = sub('\\(Intercept\\)','Intercept',colnames(design_mat))
    
# return
    result= list()
    result[['design_mat']]= design_mat

# get all possible names
    if(return_contrasts){
        possible_contrasts = sapply(sample_anno_columns, function(cova){ 
                cname = remove_na(unique(design_model_df[,cova]))
                sapply(cname, function(cn){
                        cova_lvl = paste0(cova,cn)
                        ifelse( cova_lvl %in% colnames(design_mat), 
                            cova_lvl, 
                            {rest= grep(cova, colnames(design_mat), value=T)
                                paste('-(',paste(rest, collapse='+'),')')
                            }
                        )
                    })
            }
        )
        possible_contrasts = unlist(possible_contrasts, recursive = F)
        
# set matrix
        contrast_mat <- makeContrasts( contrasts=possible_contrasts, levels=design_mat)
        colnames(contrast_mat)= gsub('\\.','',names(possible_contrasts))
        result[['contrast_mat']]= contrast_mat
    }
# return
    return(result)
}


#get_proteome_design_matrix('35791')

