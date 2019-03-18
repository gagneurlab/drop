#'---
#' title: Compare Expression across Growth Media
#' author: mumichae
#' wb:
#'  input:
#'   - counts_glu: '`sm config["PROC_RESULTS"] + "/v29_overlap/counts/fib_ss/filtered_counts.Rds"`'
#'   - counts_gal: '`sm config["PROC_RESULTS"] + "/v29_overlap/counts/gal/filtered_counts.Rds"`'
#' output: 
#'   html_document:
#'    code_folding: show
#'    code_download: TRUE
#'---

#+ echo=F
saveRDS(snakemake, 'tmp/medium.Rds')
# snakemake <- readRDS('tmp/medium.Rds')
suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(cowplot)
    library(LSD)
})

counts_glu <- as.data.table(readRDS(snakemake@input$counts_glu), keep.rownames = 'gene_name')
counts_gal <- as.data.table(readRDS(snakemake@input$counts_gal), keep.rownames = 'gene_name')

sample_anno <- fread(snakemake@config$SAMPLE_ANNOTATION)
fibr2rna <- sample_anno[, .(FIBROBLAST_ID, RNA_ID)]

common_ids <- merge(fibr2rna[RNA_ID %in% colnames(counts_gal)], 
      fibr2rna[RNA_ID %in% colnames(counts_glu)], 
      by = 'FIBROBLAST_ID', 
      suffixes = c(".gal", ".glu"))

glu_dt <- melt(counts_glu[, c('gene_name', unique(common_ids$RNA_ID.glu)), with = F], 
               id.vars = 'gene_name', variable.name = 'RNA_ID', value.name = 'count')
glu_dt <- merge(glu_dt, fibr2rna, by = 'RNA_ID')

gal_dt <- melt(counts_gal[, c('gene_name', unique(common_ids$RNA_ID.gal)), with = F], 
               id.vars = 'gene_name', variable.name = 'RNA_ID', value.name = 'count')
gal_dt <- merge(gal_dt, fibr2rna, by = 'RNA_ID')

combined_counts <- merge(glu_dt, gal_dt, by = c('gene_name', 'FIBROBLAST_ID'), suffixes = c(".gal", ".glu"))

heatscatter(combined_counts$count.glu, combined_counts$count.gal, log = 'xy')
abline(0,1)

boxplot(log(combined_counts$count.glu+1), log(combined_counts$count.gal+1), xaxt = 'n')
axis(1, 1:2, c('glu', 'gal'))
wilcox.test(log(combined_counts$count.glu+1), log(combined_counts$count.gal+1))

