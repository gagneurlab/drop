# R function
#
# Author: baderd
###############################################################################


#' get_protein_limma_differential_expression
#' 
#' Compute differential expression for log2 protein intensities 
#' with moderated t-test from limma package 
#' using customized normalization via specified design matrix.
#' 
#' @author baderda
#' 
#' @param protein_intensity_mat numeric matrix; log2 protein intensity values
#' @param design_matrix matrix; output of `get_proteome_design_matrix()` 
#'      from mitomultiomics project, which uses `model.matrix()` internally
#' @param coeff_patient character; the column of design_matrix you want to compute 
#'      fold change and pvalues for
#'
#'
get_protein_limma_differential_expression <- function(
    protein_intensity_mat,
    design_matrix,
    coeff_patient = "FIBROBLAST_IDcase"
){
    require(limma)
    require(data.table)
    
    #' initial Fit
    lfit= lmFit( protein_intensity_mat, design = design_matrix)
    
    #' * perform moderated t-Test: compute P-values with empirical Bayes
    #' 
    mod_tstats <- eBayes( contrasts.fit( lfit, coefficients = coeff_patient) )
    
    #' * get fold changes of interest
    #' 
    prot_limma_de_res= as.data.table(topTable(
        mod_tstats, coef = coeff_patient,
        number=Inf, confint=T, sort.by='none', genelist=rownames(mod_tstats)
    ))
    setnames(prot_limma_de_res, c(
        'gene', 'prot_log2fc', 'prot_ci_left', 'prot_ci_right', 'prot_baseMean',
        'prot_mod_t_stat', 'prot_pvalue', 'prot_fdr', 'prot_logodds'
    ))
    
    return(prot_limma_de_res)
}
