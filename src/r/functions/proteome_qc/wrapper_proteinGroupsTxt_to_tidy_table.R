#'---
#' title: R script
#' author: Daniel Bader
#' output: 
#'   html_document:
#'     toc_float: yes
#'     toc: yes
#'---
#'

wrapper_proteinGroupsTxt_to_tidy_table <- function(
    file_protein_groups_txt,
    intensity_column_pattern = "LFQ.intensity.",
    column_ids = c("Protein.IDs", "Gene.names"),
    column_protein_id='PROTEIN_ID',
    column_gene_id='GENE_NAME', 
    column_sample_id = "PROTEOME_ID", 
    column_intensity = "LFQ_INTENSITY",
    verbose = T, 
    replace_zero_with_NA = TRUE,
    remove_proteins_without_gene_name=TRUE,
    keep_first_id_only=TRUE
){

    # read input file to tidy table
    pdt <- read_data_table_from_proteinGroupsTxt(
        file_protein_groups_txt, 
        intensity_column_pattern = intensity_column_pattern, 
        column_ids = column_ids, 
        column_sample_id = column_sample_id, 
        column_intensity = column_intensity,
        replace_zero_with_NA = replace_zero_with_NA
    )
    
    
#' Filter protein expression by gene properties
#' 
    pdt <- filter_proteome_by_gene_properties(
        pdt, 
        column_protein_id=column_protein_id,
        column_gene_id=column_gene_id,
        column_intensity=column_intensity,
        remove_proteins_without_gene_name=remove_proteins_without_gene_name,
        keep_first_id_only=keep_first_id_only,
        verbose=verbose
    )
    
    
#' Compute NA frequencies
#' 
    compute_na_frequency(pdt, 
        column_id = column_gene_id,
        column_intensity=column_intensity
    )
    compute_na_frequency(pdt, 
        column_id = column_sample_id,
        column_intensity=column_intensity
    )
    
    # remove NA-only genes
    pdt <- pdt[NA_FREQ_BY_GENE_NAME!=1]
    
    # END
    return(pdt)
}





