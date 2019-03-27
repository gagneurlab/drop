#'---
#' title: Annotated Variants
#' author: mumichae
#' wb:
#'  input:
#'   - vcf_dts: '`sm expand(config["RAW_DATA"] + "/{vcf}/exomicout/paired-endout/processedData/vep_anno_{vcf}_uniq_dt.Rds", vcf=set(config["vcfs"]))`'
#'  output:
#'   - variant_dt:  '`sm config["PROC_RESULTS"] + "/process_vcf/variant_dt.Rds"`'
#'  threads: 20
#' output: 
#'   html_document:
#'    code_folding: hide
#'    code_download: TRUE
#'---

saveRDS(snakemake, 'tmp/variant_annotation_total.Rds')
# snakemake <-  readRDS('tmp/variant_annotation_total.Rds')

suppressPackageStartupMessages({
    library(data.table)
    library(ggpubr)
    library(dplyr)
    library(ggplot2)
    library(cowplot)
    library(BiocParallel)
    })
source("Scripts/_functions/filter_sets.R")
source('Scripts/_functions/gene_annotation/add_gene_info_cols.R')

register(MulticoreParam(snakemake@threads))
all_vcfs_list <- bplapply(snakemake@input$vcf_dts, function(f) {
    dt <- readRDS(f)
    dt <- filter_vcf_quality(dt)
    summary_dt <- data.table(
        sample = dt$sample,
        quality = T,
        exonic = dt$var_id %in% filter_exonic(dt)$var_id,
        prot_effect = dt$var_id %in% filter_prot_effect(dt)$var_id,
        rare = dt$var_id %in% filter_rare(dt)$var_id,
        compound_heterozygous = dt$var_id %in% filter_only_compound_heterozygous(dt)$var_id,
        homozygous = dt$var_id %in% filter_homozygous(dt)$var_id,
        potential_biallelic = dt$var_id %in% filter_potential_biallelic(dt)$var_id
    )
    summary_dt[, homozygous_filtered := exonic & prot_effect & rare & homozygous]
    summary_dt[, compound_heterozygous_filtered := exonic & prot_effect & rare & compound_heterozygous]
    summary_dt[, potential_biallelic_filtered := exonic & prot_effect & rare & potential_biallelic]
    summary_dt <- cbind(dt[, .(var_id, chr, pos, ref, alt, mstype, mtype, MAX_AF, gnomAD_AF, hgncid, sift1, pph1, CADD_raw, CADD_phred)], summary_dt)
    
    # add mito
    summary_dt <- add_hans_class(summary_dt, gene_name_col = 'hgncid', return_all_info = F)
    summary_dt <- add_mitocarta_col(summary_dt, gene_name_col = 'hgncid')
    
    # summary_dt[, mito := F]
    # summary_dt[MITO_DISEASE_GENE == T, mito := T]
    summary_dt[, rare_mito := F]
    summary_dt[exonic == T & prot_effect == T & rare == T & MITO_DISEASE_GENE == T, rare_mito := T]
    
    summary_dt
})

all_vcfs_dt <- rbindlist(all_vcfs_list)
saveRDS(all_vcfs_dt, snakemake@output$variant_dt)
sapply(c("missense", "synonymous", "splice", "unstop", "frame-shift", "unstart", "stop", "stop_retain", "coding"), function(var_type){
  sub_vt <- all_vcfs_dt[mstype == var_type]
  sub_vt[, c("var_id", "quality") := NULL]
  saveRDS(sub_vt, paste0('/s/project/genetic_diagnosis/processed_results/process_vcf/', var_type, '_variant_dt.Rds'))
})


filters <- c('quality', 'exonic', 'prot_effect', 'rare', 'compound_heterozygous', 'homozygous', 'potential_biallelic',
             'potential_biallelic_filtered', 'rare_mito', 'homozygous_filtered', 'compound_heterozygous_filtered')
single_filters <- c('quality', 'exonic', 'prot_effect', 'rare', 'compound_heterozygous', 'homozygous', 'potential_biallelic')
total_dt <- melt(all_vcfs_dt, measure.vars = filters, variable.name = 'filter', value.name = 'filtered')
total_dt <- total_dt[filtered == T, -'filtered']
total_dt[, single_filter := F]
total_dt[filter %in% single_filters, single_filter := T]
setkey(total_dt, chr, pos, ref, alt)

#' ## Revisit Filters
filtered_dt <- total_dt[, .(variant_count = .N), by = c('filter', 'sample', 'single_filter')]
stat_box_data <- function(y, upper_limit = max(filtered_dt$variant_count)) {
    data.frame(y = upper_limit, label = round(median(y), 1))
}
#+ filtered, fig.height=8, fig.width=10
ggplot(filtered_dt, aes(reorder(filter, -variant_count), variant_count, col = single_filter)) +
    geom_boxplot() +
    stat_summary(fun.data = stat_box_data, geom = "text", vjust = -0.5) +
    coord_trans(y="log10") +
    labs(x = "Filters", y = 'Number of variants per sample', title = 'Filtering over all annotated variants') +
    grids() +
    theme(axis.text.x = element_text(angle=45, vjust = 0.5),
          legend.position = 'bottom') +
    scale_x_discrete(breaks=levels(filtered_dt[,filter]),
                     labels=gsub("_" , "\n" , levels(filtered_dt[,filter]))) +
    scale_color_brewer(palette = 'Set1')


#' ### Get Outliers
max_per_filter <- filtered_dt[, max(variant_count), by = filter]$V1

#' Top 5 outliers
filtered_dt[, median := median(variant_count), by = filter]
setorder(filtered_dt, -variant_count)
filtered_dt[, .SD[1:5], by = filter]


#' ## AF Distributions
total_dt <- total_dt[!duplicated(total_dt, by = c('chr', 'pos', 'ref', 'alt', 'filter'))]
# mark NA
total_dt[, MAX_AF_na := MAX_AF]
total_dt[is.na(MAX_AF), MAX_AF_na := -0.1]

af_dt <- total_dt[filter != 'rare' & filter %in% single_filters]
hist(af_dt[filter == 'quality', MAX_AF_na])

ggplot(af_dt, aes(MAX_AF_na, fill = filter)) +
    geom_histogram(bins = 30) +
    facet_wrap(~filter) +
    labs(x = 'maximum allele frequency') +
    theme(legend.position = 'none') +
    scale_fill_brewer(palette = "Set2")


#' ## Mito Genes
rare_variants <- unique(total_dt[filter == 'potential_biallelic', .(var_id, MITO_DISEASE_GENE, MITOCARTA)])
rare_variants <- melt(rare_variants, id.vars = 'var_id', variable.name = 'mito_class', value.name = 'mito')

ggplot(rare_variants[, .N, by = c('mito_class', 'mito')], aes(mito, N, fill = mito)) +
    geom_bar(stat = 'identity') +
    facet_wrap(.~mito_class) +
    geom_text(aes(label=N), vjust=-0.3) +
    scale_fill_brewer(palette = "Set1") +
    scale_y_log10() +
    labs(y = 'Number of unique rare variants') +
    theme(legend.position = 'none')

