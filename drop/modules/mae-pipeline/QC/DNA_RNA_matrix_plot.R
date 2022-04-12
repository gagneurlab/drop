#'---
#' title: DNA-RNA matching matrix
#' author: vyepez
#' wb:
#'  log:
#'   - snakemake: '`sm str(tmp_dir / "MAE" / "{dataset}" / "QC_matrix_plot.Rds")`'
#'  input:
#'    - mat_qc: '`sm cfg.getProcessedResultsDir() + 
#'               "/mae/{dataset}/dna_rna_qc_matrix.Rds"`'
#'  output:
#'    - wBhtml: '`sm config["htmlOutputPath"] + "/MonoallelicExpression/QC/{dataset}.html"`'
#'  type: noindex
#'---

#+echo=F
saveRDS(snakemake, snakemake@log$snakemake)

suppressPackageStartupMessages({
  library(reshape2)
  library(data.table)
  library(ggplot2)
})


#'
#' ## Plot DNA - RNA matching matrix
qc_mat <- readRDS(snakemake@input$mat_qc)
# hist(qc_mat, xlab = '% of overlapping variants from DNA and RNA', main = '')
melt_mat <- as.data.table(reshape2::melt(qc_mat))

identityCutoff <- .85

ggplot(melt_mat, aes(value)) + geom_histogram(fill = 'cadetblue4', binwidth = 0.05, center = .025) + 
  theme_bw(base_size = 14) + 
  labs(x = 'Proportion of matching DNA-RNA variants', y = 'DNA-RNA combinations') + 
  scale_y_log10() + annotation_logticks(sides = "l") + 
  expand_limits(x=c(0,1)) +
  geom_vline(xintercept=identityCutoff, linetype='dashed', color = 'firebrick')


#' ## Identify matching samples

#' Number of samples: `r nrow(qc_mat)`
#' 
#' Number of samples that match with another: `r length(qc_mat[qc_mat > identityCutoff])`
#'
#' Median of proportion of matching variants in matching samples: `r round(median(qc_mat[qc_mat > identityCutoff]), 2)`
#'
#' Median of proportion of matching variants in not matching samples: `r round(median(qc_mat[qc_mat < identityCutoff]), 2)`
#'
#' **Considerations:**
#' On our experience, the median of the proportion of matching variants in matching samples is around 0.95,
#' and the median of the proportion of matching variants in not matching samples is around 0.58.
#' Sometimes we do see some values between 0.7 - 0.85. That could mean that the DNA-RNA combination is 
#' not from the same person, but from a relative. It could also be due to a technical error. For those cases, 
#' check the following:
#' 
#' * RNA sequencing depth (low seq depth that can lead to variants not to be found in the RNA)
#' * Number of variants (too many variants called due to sequencing errors)
#' * Ratio of heterozygous/homozygous variants (usually too many called variants means too many heterozygous ones)
#' * Is the sample a relative of the other?
#' 

sa <- fread(snakemake@config$sampleAnnotation)[, .(DNA_ID, RNA_ID)]
sa[, ANNOTATED_MATCH := TRUE]
colnames(melt_mat)[1:2] <- c('DNA_ID', 'RNA_ID')

#' ### Samples that were annotated to match but do not 
false_matches <- merge(sa, melt_mat, by = c('DNA_ID', 'RNA_ID'), 
                       sort = FALSE, all.x = TRUE)
DT::datatable(false_matches[value < identityCutoff])

#' ### Samples that were not annotated to match but actually do
false_mismatches <- merge(melt_mat, sa, by = c('DNA_ID', 'RNA_ID'), 
                          sort = FALSE, all.x = TRUE)
DT::datatable(false_mismatches[is.na(ANNOTATED_MATCH) & value > identityCutoff])

