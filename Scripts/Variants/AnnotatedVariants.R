#'---
#' title: Annotated Variants
#' author: mumichae
#' wb:
#'  input:
#'   - vcf_dts: '`sm expand(config["RAW_DATA"] + "/{vcf}/exomicout/paired-endout/processedData/vep_anno_{vcf}_uniq_dt.Rds", vcf=set(config["vcfs"]))`'
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
        compound_heterzygous = dt$var_id %in% filter_only_compound_heterzygous(dt)$var_id,
        homozygous = dt$var_id %in% filter_homozygous(dt)$var_id,
        potential_biallelic = dt$var_id %in% filter_potential_biallelic(dt)$var_id
    )
    summary_dt[, all_filters := exonic & prot_effect & rare & compound_heterzygous 
               & homozygous & potential_biallelic]
    cbind(dt[, .(var_id, chr, pos, ref, alt, MAX_AF, gnomAD_AF, hgncid, sift1, pph1)], summary_dt)
})

filters <- c('quality', 'exonic', 'prot_effect', 'rare', 'compound_heterzygous', 'homozygous', 'potential_biallelic', 'all_filters')
total_dt <- melt(rbindlist(all_vcfs_list), measure.vars = filters, variable.name = 'filter', value.name = 'filtered')
total_dt <- total_dt[filtered == T, -'filtered']
setkey(total_dt, chr, pos, ref, alt)

#' ## Revisit Filters
filtered_dt <- total_dt[, .(variant_count = .N), by = c('filter', 'sample')]
ggplot(filtered_dt,
       aes(reorder(filter, -variant_count), variant_count)) +
    geom_boxplot() +
    grids() +
    labs(x = "Filters", y = 'Number of variants per sample') +
    scale_y_log10() + 
    # theme(axis.text.x = element_text(angle=45)) + 
    scale_x_discrete(breaks=levels(filtered_dt[,filter]),
                     labels=gsub("_" , "\n" , levels(filtered_dt[,filter])))


#' Get Outliers
max_per_filter <- filtered_dt[, max(variant_count), by = filter]$V1
filtered_dt[variant_count %in% max_per_filter]

# top 5
filtered_dt[, median := median(variant_count), by = filter]
setorder(filtered_dt, -variant_count)
filtered_dt[, .SD[1:5], by = filter]


#' ## AF Distributions
# mark NA
total_dt[, MAX_AF_na := MAX_AF]
total_dt[is.na(MAX_AF), MAX_AF_na := -0.1]
# keep only unique and omit MAF filters
total_dt <- total_dt[!duplicated(total_dt, by = c('chr', 'pos', 'ref', 'alt', 'filter'))]
af_dt <- total_dt[filter != 'rare' & filter != 'all_filters']

hist(af_dt$MAX_AF_na)

ggplot(af_dt, aes(MAX_AF_na, fill = filter)) +
    geom_histogram() +
    facet_wrap(~filter) +
    scale_fill_brewer(palette = "Set2")

ggplot(af_dt, aes(MAX_AF_na, stat(count), col = filter)) +
    geom_density() +
    scale_color_brewer(palette = "Set2")


#' ## Mito Genes
source('Scripts/_functions/gene_annotation/add_gene_info_cols.R')
gene_info_dt <- add_all_gene_info(total_dt, gene_name_col = 'hgncid', dis_genes = F)
# gene_info_dt <- melt(gene_info_dt, measure.vars = c('HANS_CLASS', 'MITOCARTA'), variable.name = 'mito', value.name = 'mito_desc')
# gene_info_dt <- na.omit(gene_info_dt, 'mito_desc')

ggplot(gene_info_dt[, .N, by = HANS_CLASS], aes(HANS_CLASS, N, fill = HANS_CLASS)) +
    geom_bar(stat = 'identity') +
    geom_text(aes(label=N), vjust=-0.3) +
    scale_y_log10()

ggplot(gene_info_dt[, .N, by = MITOCARTA], aes(MITOCARTA, N, fill = MITOCARTA)) +
    geom_bar(stat = 'identity') +
    geom_text(aes(label=N), vjust=-0.3) +
    scale_fill_brewer(palette = 'Set1') +
    scale_y_log10()

