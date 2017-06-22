#'---
#' title: R script
#' author: Daniel Bader
#' output: 
#'   html_document:
#'     toc_float: yes
#'     toc: yes
#'---
#'




compute_na_frequency <- function(
    proteome_data_table,
    column_id,
    new_column_label=paste0("NA_FREQ_BY_",column_id),
    column_intensity='LFQ_INTENSITY'
){
    proteome_data_table[,
        eval(new_column_label) := sum(is.na(get(column_intensity)))/.N, 
        by=get(column_id)
    ]
}






