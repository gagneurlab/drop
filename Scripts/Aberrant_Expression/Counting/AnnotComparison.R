#'---
#' title: Compare Annotations v19 and v29
#' author: Michaela Muller, vyepez
#' wb:
#'  input: 
#'   - filtered_v19: '/s/project/genetic_diagnosis/processed_results/v19/outrider/ods_unfitted.Rds'
#'   - filtered_v29: '/s/project/genetic_diagnosis/processed_results/v29/outrider/ods_unfitted.Rds'
#' output: 
#'   html_document:
#'    code_folding: show
#'    code_download: TRUE
#'---

#+ echo=F
saveRDS(snakemake, "tmp/count_analysis.snakemake")
# snakemake <- readRDS("tmp/count_analysis.snakemake")
suppressPackageStartupMessages({
    library(SummarizedExperiment)
    library(data.table)
    library(tidyr)
    library(ggplot2)
    library(magrittr)
})
source("Scripts/_functions/gene_annotation/add_gene_info_cols.R")

#' ## Goal
#' Compare both v19 and v29 annotations by computing the total number of genes that
#' exist in both gtf files, the number of genes with at least 1 count in 1 sample,
#' the number of genes that passed the OUTRIDER filter.
#' Also consider protein coding and mito disease genes.

#' ## Read and combine the rowData from both annotations
filtv19 <- readRDS(snakemake@input$filtered_v19)
rd19 <- rowData(filtv19) %>% as.data.table()
rd19[, version := 'v19']

filtv29 <- readRDS(snakemake@input$filtered_v29)
rd29 <- rowData(filtv29) %>% as.data.table()
rd29[, version := 'v29']
rd <- rbind(rd19, rd29)

# Add mito disease genes information
rd[, gene_name_unique := toupper(gene_name_unique)]
rd <- add_hans_class(rd, gene_name_col = "gene_name_unique", return_all_info = FALSE)


# List of special genes that Robert asked to check
genes_to_check <- c("NDUFAF5", "NDUFAF6",
                    "NDUFA2", "NDUFA7", "NDUFA3", "NDUFA11", "NDUFA13",
                    "DNAJC30", "GCDH", "MOCS1",
                    "MRPL12", "MRPL30", "MRPL38", "MRPL45",
                    "MRPS17", "MRPS21",
                    "MSRB3", "MTG1", "RARS2", "SARS2",
                    "SDHAF2", "SLC25A10", "SLC25A26",
                    "TIMM9", "TIMM10B", "TIMM13", "TIMM23",
                    "TK2", "TOMM5", "TSFM", "ACACA", "ACAD11", "FAHD1", "GATC")
length(genes_to_check)

rd[gene_name_unique %in% genes_to_check, gene_to_check := TRUE]

#' ## Create comparison table
comp_dt <- data.table(annotation = c("v19", "v29"))
# How many genes are in the GTF file?
comp_dt[, total_all := nrow(rd[version == annotation]), by = 1:nrow(comp_dt)]
comp_dt[, total_pc := nrow(rd[version == annotation & gene_type == 'protein_coding']), by = 1:nrow(comp_dt)]
comp_dt[, total_mito := nrow(rd[version == annotation & MITO_DISEASE_GENE == TRUE]), by = 1:nrow(comp_dt)]
comp_dt[, total_sp := nrow(rd[version == annotation & gene_to_check == TRUE]), by = 1:nrow(comp_dt)]

# How many genes have at least 1 read for all samples?
comp_dt[, counted_all := sum(rd[version == annotation, counted1sample]), by = 1:nrow(comp_dt)]
comp_dt[, counted_pc := sum(rd[version == annotation & gene_type == 'protein_coding', counted1sample]), by = 1:nrow(comp_dt)]
comp_dt[, counted_mito := sum(rd[version == annotation & MITO_DISEASE_GENE == TRUE, counted1sample]), by = 1:nrow(comp_dt)]
comp_dt[, counted_sp := sum(rd[version == annotation & gene_to_check == TRUE, counted1sample]), by = 1:nrow(comp_dt)]

# How many genes pass the OUTRIDER filter?
comp_dt[, passedFilter_all := sum(rd[version == annotation, passedFilter]), by = 1:nrow(comp_dt)]
comp_dt[, passedFilter_pc := sum(rd[version == annotation & gene_type == 'protein_coding', passedFilter]), by = 1:nrow(comp_dt)]
comp_dt[, passedFilter_mito := sum(rd[version == annotation & MITO_DISEASE_GENE == TRUE, passedFilter]), by = 1:nrow(comp_dt)]
comp_dt[, passedFilter_sp := sum(rd[version == annotation & gene_to_check == TRUE, passedFilter]), by = 1:nrow(comp_dt)]
comp_dt

#' ## Tidy up data and plotting
mt <- melt(comp_dt)
mt <- separate(mt, col = 'variable', into = c('group', 'category'), sep = "_")
mt[, prop := value/max(value), by = category]
mt[, group := factor(group, levels = c("total", "counted", "passedFilter"))]
mt[, annotation := factor(annotation, levels = c("v29", "v19"))]
mt[category == 'mito', category := 'Mito_disease']
mt[category == 'pc', category := 'protein_coding']

library(ggthemes)
ggplot(mt[category == 'all'], aes(group, prop, fill = annotation)) + geom_bar(stat = 'identity', position = 'dodge') + 
    geom_text(aes(label = value),  position = position_dodge(width = .8), vjust = -.5) + 
    theme_bw(base_size = 14) + scale_fill_brewer(palette="BuPu")

#+ fig.width=9, fig.height=10
ggplot(mt, aes(group, prop, fill = annotation)) + geom_bar(stat = 'identity', position = 'dodge') + 
    geom_text(aes(label = value),  position = position_dodge(width = .8), vjust = -.3) + 
    theme_bw(base_size = 14) + scale_fill_brewer(palette="BuPu") + facet_wrap(~category)

#'
#' Mito disease genes with no counts:
rd[version == 'v29' & MITO_DISEASE_GENE == T & counted1sample == F, gene_name_unique] %>% sort

#' Mito disease genes that didn't pass the OUTRIDER filter:
rd[version == 'v29' & MITO_DISEASE_GENE == T & counted1sample == T & passedFilter == F, gene_name_unique] %>% sort

#' Robert's special genes with no counts:
rd[version == 'v29' & gene_to_check == T & counted1sample == F, gene_name_unique] %>% sort

#' Robert's special genes that didn't pass the OUTRIDER filter:
rd[version == 'v29' & gene_to_check == T & counted1sample == T & passedFilter == F, gene_name_unique] %>% sort



