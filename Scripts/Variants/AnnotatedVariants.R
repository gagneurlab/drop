#'---
#' title: Annotated Variants
#' author: mumichae
#' wb:
#'  input:
#'   - vcf_dts: '`sm expand(config["RAW_DATA"] + "/{vcf}/exomicout/paired-endout/processedData/vep_anno_{vcf}_uniq_dt.Rds", vcf=config["vcfs"])`'
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
    library(ggplot2)
    library(cowplot)
    library(BiocParallel)
    })
source("Scripts/_functions/filter_sets.R")

register(MulticoreParam(snakemake@threads))
all_vcfs_list <- bplapply(snakemake@input$vcf_dts, function(f) {
    dt <- readRDS(f)
    dt <- filter_vcf_quality(dt)
    dt[, id := 1:.N]
    data.table(
        sample = dt$sample,
        quality = T,
        exonic = dt$id %in% filter_exonic(dt)$id,
        prot_effect = dt$id %in% filter_prot_effect(dt)$id,
        rare = dt$id %in% filter_rare(dt)$id,
        compound_heterzygous = dt$id %in% filter_only_compound_heterzygous(dt)$id,
        homozygous = dt$id %in% filter_homozygous(dt)$id,
        potential_biallelic = dt$id %in% filter_potential_biallelic(dt)$id,
        all_filters = dt$id %in% 
            filter_vcf_quality(
            filter_exonic(
            filter_prot_effect(
            filter_rare(
            filter_only_compound_heterzygous(
            filter_homozygous(
            filter_potential_biallelic(dt)))))))$id,
        MAX_AF = dt$max_maf #MAX_AF # need to recompute the data.tables to update MAF, currently "max_maf"
    )
})

filters <- c('quality', 'exonic', 'prot_effect', 'rare', 'compound_heterzygous', 'homozygous', 'potential_biallelic')
total_dt <- melt(rbindlist(all_vcfs_list), measure.vars = filters, variable.name = 'filter', value.name = 'filtered')

#' ## Revisit Filters
filter_dt <- total_dt[filtered == T, .(variant_count = .N), by = c('filter', 'sample')]
ggplot(filter_dt,
       aes(reorder(filter, variant_count), variant_count)) +
    geom_boxplot() +
    labs(x = "Filters", y = 'Number of variants per sample') +
    coord_flip()

#' ## MAF Distributions
ggplot(total_dt[filtered == T], aes(reorder(filter, MAX_AF), MAX_AF)) +
    geom_boxplot()

ggplot(total_dt[filtered == T], aes(MAX_AF)) +
    geom_histogram() +
    facet_wrap(~filter)

ggplot(total_dt[filtered == T], aes(MAX_AF, col = filter)) +
    geom_density() +
    scale_y_log10()
