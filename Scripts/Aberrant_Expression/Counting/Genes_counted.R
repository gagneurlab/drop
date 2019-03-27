#'---
#' title: Genes counted in different tissues
#' author: vyepez
#' wb:
#'  input: 
#'   - fib_ss: '/s/project/genetic_diagnosis/processed_results/v29_overlap/outrider/fib_ss/ods_unfitted.Rds'
#'   - blood: '/s/project/genetic_diagnosis/processed_results/v29_overlap/outrider/blood/ods_unfitted.Rds'
#'   - add_gene_info_script: 'Scripts/_functions/gene_annotation/add_gene_info_cols.R'
#' output: 
#'   html_document:
#'    code_folding: show
#'    code_download: TRUE
#'---

#+ echo=F
saveRDS(snakemake, "tmp/genes_tissues.snakemake")
# snakemake <- readRDS("tmp/genes_tissues.snakemake")
suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(data.table)
  library(tidyr)
  library(ggplot2)
  library(magrittr)
})
source(snakemake@input$add_gene_info_script)


#' the number of genes that passed the OUTRIDER filter.
#' Also consider protein coding and mito disease genes.
#'
#' ## Read and combine the rowData from both annotations
ods_ss <- readRDS(snakemake@input$fib_ss)
row.names(ods_ss) <- rowData(ods_ss)$gene_name_unique
rd_ss <- rowData(ods_ss) %>% as.data.table()
rd_ss[, group := 'fib_ss']

ods_blood <- readRDS(snakemake@input$blood)
row.names(ods_blood) <- rowData(ods_blood)$gene_name_unique
rd_blood <- rowData(ods_blood) %>% as.data.table()
rd_blood[, group := 'blood']

rd <- rbind(rd_ss, rd_blood, fill = T)

# Add mito disease genes information
rd[, gene_name_unique := toupper(gene_name_unique)]
rd <- add_hans_class(rd, gene_name_col = "gene_name_unique", return_all_info = FALSE)
rd <- add_omim_cols(rd, gene_name_col = "gene_name_unique", return_all_info = FALSE)

#' ## Create comparison table
comp_dt <- data.table(tissue = unique(rd$group))
# How many genes are in the GTF file?
comp_dt[, total_all := nrow(rd[group == tissue]), by = 1:nrow(comp_dt)]
comp_dt[, total_pc := nrow(rd[group == tissue & gene_type == 'protein_coding']), by = 1:nrow(comp_dt)]
comp_dt[, total_mito := nrow(rd[group == tissue & MITO_DISEASE_GENE == TRUE]), by = 1:nrow(comp_dt)]
comp_dt[, total_omim := nrow(rd[group == tissue & OMIM_gene == TRUE]), by = 1:nrow(comp_dt)]

# How many genes have at least 1 read for all samples?
comp_dt[, counted_all := sum(rd[group == tissue, counted1sample]), by = 1:nrow(comp_dt)]
comp_dt[, counted_pc := sum(rd[group == tissue & gene_type == 'protein_coding', counted1sample]), by = 1:nrow(comp_dt)]
comp_dt[, counted_mito := sum(rd[group == tissue & MITO_DISEASE_GENE == TRUE, counted1sample]), by = 1:nrow(comp_dt)]
comp_dt[, counted_omim := sum(rd[group == tissue & OMIM_gene == TRUE, counted1sample]), by = 1:nrow(comp_dt)]

# How many genes pass the OUTRIDER filter?
comp_dt[, passedFilter_all := sum(rd[group == tissue, passedFilter]), by = 1:nrow(comp_dt)]
comp_dt[, passedFilter_pc := sum(rd[group == tissue & gene_type == 'protein_coding', passedFilter]), by = 1:nrow(comp_dt)]
comp_dt[, passedFilter_mito := sum(rd[group == tissue & MITO_DISEASE_GENE == TRUE, passedFilter]), by = 1:nrow(comp_dt)]
comp_dt[, passedFilter_omim := sum(rd[group == tissue & OMIM_gene == TRUE, passedFilter]), by = 1:nrow(comp_dt)]

DT::datatable(comp_dt, style = 'bootstrap')

#' ## Tidy up data and plotting
mt <- melt(comp_dt)
mt <- separate(mt, col = 'variable', into = c('group', 'category'), sep = "_")
mt[, prop := value/max(value, na.rm =T), by = category]
mt[, group := factor(group, levels = c("total", "counted", "passedFilter"))]
mt[, tissue := factor(tissue, levels = c("blood", "fib_ss"))]
mt[category == 'mito', category := 'Mito_disease']
mt[category == 'pc', category := 'protein_coding']
mt <- na.omit(mt)


#+ fig.width=13, fig.height=12
ggplot(mt, aes(group, prop, fill = tissue)) + geom_bar(stat = 'identity', position = 'dodge') + 
  geom_text(aes(label = value),  position = position_dodge(width = .8), vjust = -.1) + 
  theme_bw(base_size = 14) + scale_fill_brewer(palette = "Set2") + facet_wrap(~category, ncol = 2)


#' Mito disease genes with no counts:
rd[MITO_DISEASE_GENE == T & counted1sample == F, .(group, gene_name_unique)]

#' Mito disease genes with counts that didn't pass the OUTRIDER filter:
rd[MITO_DISEASE_GENE == T & counted1sample == T & passedFilter == F, .(group, gene_name_unique)] 

#' Number of nuclear mito disease genes that passed the OUTRIDER filter:
mito_nuclear_filt <- rd[MITO_DISEASE_GENE == T & gene_type == 'protein_coding' & passedFilter == T, .N, by = group]
mito_nuclear_filt[, Nprop := round(N/rd[MITO_DISEASE_GENE == T & gene_type == 'protein_coding', .N]*uniqueN(rd$group), 3)]
mito_nuclear_filt

#' Genes that pass the filter in fibroblasts and not in blood
setdiff(rd[MITO_DISEASE_GENE == T & passedFilter == T & group == 'fib_ss', gene_name_unique], 
        rd[MITO_DISEASE_GENE == T & passedFilter == T & group == 'blood', gene_name_unique])

#' Genes that pass the filter in blood and not in fibroblasts
setdiff(rd[MITO_DISEASE_GENE == T & passedFilter == T & group == 'blood', gene_name_unique], 
        rd[MITO_DISEASE_GENE == T & passedFilter == T & group == 'fib_ss', gene_name_unique])

