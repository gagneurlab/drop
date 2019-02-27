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
summary_list <- bplapply(snakemake@input$vcf_dts, function(f) {
    dt <- readRDS(f)
    dt <- filter_vcf_quality(dt)
    data.table(
        id = dt$sample[[1]],
        quality = nrow(dt),
        exome = nrow(filter_exome(dt)),
        prot_effect = nrow(filter_prot_effect(dt)),
        rare = nrow(filter_rare(dt)),
        compound_heterzygous = nrow(filter_only_compound_heterzygous(dt)),
        homozygous = nrow(filter_homozygous(dt)),
        potential_biallelic = nrow(filter_potential_biallelic(dt)),
        
        all_filters = nrow(
            filter_vcf_quality(
            filter_exome(
            filter_prot_effect(
            filter_rare(
            filter_only_compound_heterzygous(
            filter_homozygous(
            filter_potential_biallelic(dt))))))))
    )
})

summary_dt <- melt(rbindlist(summary_list), id.vars = 'id',  variable.name = 'filter', value.name = 'variant_count')

ggplot(summary_dt, aes(reorder(filter, variant_count), variant_count)) +
    geom_boxplot() +
    labs(x = "Filters", y = 'Number of variants') +
    coord_flip()

