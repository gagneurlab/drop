#'---
#' title: Compare Annotations v19 and v29
#' author: Michaela Muller, vyepez
#' wb:
#'  input: 
#'   - filtered_v19: '/s/project/genetic_diagnosis/processed_results/v19/outrider/ods_unfitted.Rds'
#'   - filtered_v29: '/s/project/genetic_diagnosis/processed_results/v29/outrider/ods_unfitted.Rds'
#'   - filtered_v29_ov: '/s/project/genetic_diagnosis/processed_results/v29_overlap/outrider/fib_ss/ods_unfitted.Rds'
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
#'
#' ## Read and combine the rowData from both annotations
filtv19 <- readRDS(snakemake@input$filtered_v19)
row.names(filtv19) <- rowData(filtv19)$gene_name_unique
rd19 <- rowData(filtv19) %>% as.data.table()
rd19[, version := 'v19']

ov29 <- readRDS(snakemake@input$filtered_v29_ov)
row.names(ov29) <- rowData(ov29)$gene_name_unique
rd29ov <- rowData(ov29) %>% as.data.table()
rd29ov[, version := 'v29ov']


filtv29 <- readRDS(snakemake@input$filtered_v29)
row.names(filtv29) <- rowData(filtv29)$gene_name_unique
rd29 <- rowData(filtv29) %>% as.data.table()
rd29[, version := 'v29']

rd <- rbind(rd19, rd29, rd29ov, fill = T)

# Add mito disease genes information
rd[, gene_name_unique := toupper(gene_name_unique)]
rd <- add_hans_class(rd, gene_name_col = "gene_name_unique", return_all_info = FALSE)
rd <- add_omim_cols(rd, gene_name_col = "gene_name_unique", return_all_info = FALSE)

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
comp_dt <- data.table(annotation = c("v19", "v29", "v29ov"))
# How many genes are in the GTF file?
comp_dt[, total_all := nrow(rd[version == annotation]), by = 1:nrow(comp_dt)]
comp_dt[, total_pc := nrow(rd[version == annotation & gene_type == 'protein_coding']), by = 1:nrow(comp_dt)]
comp_dt[, total_mito := nrow(rd[version == annotation & MITO_DISEASE_GENE == TRUE]), by = 1:nrow(comp_dt)]
comp_dt[, total_special := nrow(rd[version == annotation & gene_to_check == TRUE]), by = 1:nrow(comp_dt)]
comp_dt[, total_omim := nrow(rd[version == annotation & OMIM_gene == TRUE]), by = 1:nrow(comp_dt)]

# How many genes have at least 1 read for all samples?
comp_dt[, counted_all := sum(rd[version == annotation, counted1sample]), by = 1:nrow(comp_dt)]
comp_dt[, counted_pc := sum(rd[version == annotation & gene_type == 'protein_coding', counted1sample]), by = 1:nrow(comp_dt)]
comp_dt[, counted_mito := sum(rd[version == annotation & MITO_DISEASE_GENE == TRUE, counted1sample]), by = 1:nrow(comp_dt)]
comp_dt[, counted_special := sum(rd[version == annotation & gene_to_check == TRUE, counted1sample]), by = 1:nrow(comp_dt)]
comp_dt[, counted_omim := sum(rd[version == annotation & OMIM_gene == TRUE, counted1sample]), by = 1:nrow(comp_dt)]
comp_dt[annotation == 'exon', c('total_all', 'total_pc', 'total_mito', 'total_special', 'total_omim', 'counted_all', 'counted_pc', 'counted_mito', 'counted_special', 'counted_omim') := NA]

# How many genes pass the OUTRIDER filter?
comp_dt[, passedFilter_all := sum(rd[version == annotation, passedFilter]), by = 1:nrow(comp_dt)]
comp_dt[, passedFilter_pc := sum(rd[version == annotation & gene_type == 'protein_coding', passedFilter]), by = 1:nrow(comp_dt)]
comp_dt[, passedFilter_mito := sum(rd[version == annotation & MITO_DISEASE_GENE == TRUE, passedFilter]), by = 1:nrow(comp_dt)]
comp_dt[, passedFilter_special := sum(rd[version == annotation & gene_to_check == TRUE, passedFilter]), by = 1:nrow(comp_dt)]
comp_dt[, passedFilter_omim := sum(rd[version == annotation & OMIM_gene == TRUE, passedFilter]), by = 1:nrow(comp_dt)]

DT::datatable(comp_dt, style = 'bootstrap')

#' ## Tidy up data and plotting
mt <- melt(comp_dt)
mt <- separate(mt, col = 'variable', into = c('group', 'category'), sep = "_")
mt[, prop := value/max(value, na.rm =T), by = category]
mt[, group := factor(group, levels = c("total", "counted", "passedFilter"))]
mt[, annotation := factor(annotation, levels = c("exon", "v19", "v29", "v29ov"))]
mt[category == 'mito', category := 'Mito_disease']
mt[category == 'pc', category := 'protein_coding']
mt <- na.omit(mt)
mt <- mt[!(annotation == 'v29ov' & group == 'total')]

library(ggthemes)
ggplot(mt[category == 'all'], aes(group, prop, fill = annotation)) + geom_bar(stat = 'identity', position = 'dodge') + 
    geom_text(aes(label = value),  position = position_dodge(width = .8), vjust = -.5) + 
    theme_bw(base_size = 14) + scale_fill_brewer(palette = "Set2")

#+ fig.width=13, fig.height=12
ggplot(mt, aes(group, prop, fill = annotation)) + geom_bar(stat = 'identity', position = 'dodge') + 
    geom_text(aes(label = value),  position = position_dodge(width = .8), vjust = -.1) + 
    theme_bw(base_size = 14) + scale_fill_brewer(palette = "Set2") + facet_wrap(~category, ncol = 2)

#'
#' Special genes detected at first
rd[version == 'v19' & gene_to_check == T, gene_name_unique] %>% sort

#' Mito disease genes with no counts:
rd[version == 'v29' & MITO_DISEASE_GENE == T & counted1sample == F, gene_name_unique] %>% sort

#' Mito disease genes with counts that didn't pass the OUTRIDER filter:
rd[version == 'v29' & MITO_DISEASE_GENE == T & counted1sample == T & passedFilter == F, gene_name_unique] %>% sort

#' Number of nuclear mito disease genes that passed the OUTRIDER filter:
rd[version == 'v29' & MITO_DISEASE_GENE == T & gene_type == 'protein_coding' & passedFilter ==T, .N]
(rd[version == 'v29' & MITO_DISEASE_GENE == T & gene_type == 'protein_coding' & passedFilter ==T, .N] / rd[version == 'v29' & MITO_DISEASE_GENE == T & gene_type == 'protein_coding', .N]) %>% round(3)

#' Number of mtDNA mito disease genes that passed the OUTRIDER filter:
rd[version == 'v29' & MITO_DISEASE_GENE == T & gene_type != 'protein_coding' & passedFilter ==T, .N]
(rd[version == 'v29' & MITO_DISEASE_GENE == T & gene_type != 'protein_coding' & passedFilter ==T, .N] / rd[version == 'v29' & MITO_DISEASE_GENE == T & gene_type != 'protein_coding', .N]) %>% round(3)
rd[version == 'v29' & MITO_DISEASE_GENE == T & gene_type != 'protein_coding' & passedFilter ==T, gene_name_unique]

#' Robert's special genes appearing on v29ov, but not on v29
rd[version == 'v29' & gene_to_check == T & counted1sample == T & passedFilter == F, gene_name_unique] %>% sort

#' Mito disease genes appearing on v29ov, but not on v29
setdiff(rd[version == 'v29ov' & MITO_DISEASE_GENE == T & passedFilter == T, gene_name_unique], 
        rd[version == 'v29' & MITO_DISEASE_GENE == T & passedFilter == T, gene_name_unique])

fp <- fpkm(filtv29)
quantile(fp["ACAD11", ], prob = .95)


            