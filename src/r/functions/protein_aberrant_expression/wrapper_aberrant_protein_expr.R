#'---
#' title: R script
#' author: Daniel Bader
#' output: 
#'   html_document:
#'     toc_float: yes
#'     toc: yes
#'---
#'



wrapper_aberrant_protein_expr <- function(
    prot_intensity,
    sample_annotation = SAMPLE_ANNOTATION,
    coln_sample_id = 'FIBROBLAST_ID',
    coln_normalization_coeff = NULL,
    p_adjust_method = 'hochberg'
){
    library(data.table)
    res <- data.table()
    
    # normalize input by size factors
    prot_log2_intensity = normalize_expression_matrix(
        prot_intensity, sizefactor = T, rowcenter = F, log2scale = T, nonzero = 0
    )
    
    
    all_sample_ids <- sapply(
        colnames(prot_log2_intensity), get_fibro_for_proteome, sample_annotation
    )
    
    for(FIB in all_sample_ids){
        mydesign= get_proteome_design_matrix(
            FIB, 
            proteome_ids = colnames(prot_log2_intensity), 
            sample_anno_columns = c(coln_sample_id, coln_normalization_coeff),
            sample_anno_dt = sample_annotation,
            binarize_fibro_id= TRUE,
            return_contrasts = FALSE
        )
        
        # check
        sample_diff <- setdiff(
            colnames(prot_log2_intensity), rownames(mydesign$design_mat)
        )
        stopifnot(length(sample_diff)==0)
        
        
        single_proteome_limma_res= get_protein_limma_differential_expression(
            prot_log2_intensity, 
            design_matrix = mydesign[['design_mat']],
            coeff_patient = paste0(coln_sample_id, 'case')
        )
        
        # hochberg_padj
        single_proteome_limma_res[,
            prot_padj := p.adjust(prot_pvalue, method = p_adjust_method)
            ]
        # add sample
        single_proteome_limma_res[,eval(coln_sample_id):= FIB]
        
        # update res
        setnames(
            single_proteome_limma_res, toupper(names(single_proteome_limma_res))
        )
        res <- rbind(res, single_proteome_limma_res)
    }
    
    
    # 
    # COMPUTE ZSCORE
    # 
    
    #' * get tidy normalized expression
    prot_log2_intensity_by_fibro <- summarize_ematrix_by_fibroblast(
        prot_log2_intensity, 'get_fibro_for_proteome', sample_annotation
    )
    prot_log2_intensity_by_fibro <- melt(
        as.data.table(prot_log2_intensity_by_fibro, keep.rownames = T), 
        id.vars='rn',
        measure.vars = colnames(prot_log2_intensity_by_fibro),
        value.name = 'NORM_LOG2_LFQ',
        variable.name = coln_sample_id
    )
    setnames(
        prot_log2_intensity_by_fibro, 
        sub('rn', 'GENE_NAME', names(prot_log2_intensity_by_fibro))
    )
    
    #' * add std deviation per gene
    prot_log2_intensity_by_fibro[, 
        SD_NORM_LOG2_LFQ := sd(NORM_LOG2_LFQ, na.rm=T), 
        by=GENE_NAME
    ]
    
    #' * merge
    #' 
    prot_aberexp <- merge(res, prot_log2_intensity_by_fibro)
    
    #' * compute Z-score
    #' 
    prot_aberexp[, PROT_ZSCORE := PROT_LOG2FC/SD_NORM_LOG2_LFQ]
    
    
    #
    return(prot_aberexp)
}





