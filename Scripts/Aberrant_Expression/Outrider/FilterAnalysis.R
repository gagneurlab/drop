#'---
#' title: Analyze filters
#' author: vyepez
#' wb:
#'  input:
#'   - ods_ss: '`sm config["PROC_RESULTS"] + "/v29_overlap/outrider/fib_ss/ods_unfitted.Rds"`'
#'   - ods_ns: '`sm config["PROC_RESULTS"] + "/v29_overlap/outrider/fib_ns/ods_unfitted.Rds"`'
#'  output:
#' output:
#'   html_document:
#'    code_folding: show
#'    code_download: TRUE
#'---

#+ echo=F
saveRDS(snakemake, "tmp/filter_analysis.snakemake")
# snakemake <- readRDS("tmp/filter_analysis.snakemake")

suppressPackageStartupMessages({
    library(OUTRIDER)
    library(ggplot2)
    library(data.table)
    library(dplyr)
})

ods_ss <- readRDS(snakemake@input$ods_ss)
ods_ns <- readRDS(snakemake@input$ods_ns)

ods_ss <- ods_ss[mcols(ods_ss)$counted1sample, ]
ods_ns <- ods_ns[mcols(ods_ns)$counted1sample, ] 

#' ## Plot filtered samples
plotFPKM(ods_ss) + theme_bw(base_size = 14) + ggtitle('Strand Specific')

plotFPKM(ods_ns) + theme_bw(base_size = 14) + ggtitle('Non Strand Specific')


#' ## Genes with FPKM > 1
# Compute FPKM and test if it's greater than 1
fpkm_ss <- fpkm(ods_ss) > 1
fpkm_ns <- fpkm(ods_ns) > 1

#' Number of genes with FPKM > 1, per sample
fpkm_ss %>% colSums2 %>% mean
fpkm_ns %>% colSums2 %>% mean

#+ echo=F
fpkm_ss %>% colSums2 %>% hist(breaks = 15, main = 'Genes with FPKM > 1 per sample \n Strand Specific')
fpkm_ns %>% colSums2 %>% hist(breaks = 15, main = 'Genes with FPKM > 1 per sample \n Non Strand Specific')


#' ## Protein coding genes with FPKM > 1
# Subset for protein coding genes
ods_ss <- ods_ss[mcols(ods_ss)$gene_type == 'protein_coding', ]
ods_ns <- ods_ns[mcols(ods_ns)$gene_type == 'protein_coding', ]

# Compute FPKM and test if it's greater than 1
fpkm_ss <- fpkm(ods_ss) > 1
fpkm_ns <- fpkm(ods_ns) > 1

#' Number of protein coding genes with FPKM > 1, per sample
fpkm_ss %>% colSums2 %>% mean
fpkm_ns %>% colSums2 %>% mean

#+ echo=F
fpkm_ss %>% colSums2 %>% hist(breaks = 15, main = 'Protein coding genes with FPKM > 1 per sample \n Strand Specific')
fpkm_ns %>% colSums2 %>% hist(breaks = 15, main = 'Protein coding genes with FPKM > 1 per sample \n Non Strand Specific')


