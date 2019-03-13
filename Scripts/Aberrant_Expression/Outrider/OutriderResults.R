#'---
#' title: Analyze OUTRIDER object and results
#' author: Michaela Mueller, vyepez
#' wb:
#'  input:
#'   - ods_ss: '`sm expand(config["PROC_RESULTS"] + "/{annotation}/outrider/fib_ss/ods.Rds", annotation=config["ANNOTATIONS"])`'
#'   - ods_ns: '`sm expand(config["PROC_RESULTS"] + "/{annotation}/outrider/fib_ns/ods.Rds", annotation=config["ANNOTATIONS"])`'
#'  output:
#'   - results: '`sm config["PROC_RESULTS"] + "/v29_overlap/outrider/OUTRIDER_results.tsv"`'
#' output:
#'   html_document:
#'    code_folding: show
#'    code_download: TRUE
#'---

#+ echo=F
saveRDS(snakemake, "tmp/outrider_analysis.snakemake")
# snakemake <- readRDS("tmp/outrider_analysis.snakemake")

suppressPackageStartupMessages({
    library(OUTRIDER)
    library(SummarizedExperiment)
    library(ggplot2)
    library(data.table)
    library(dplyr)
})

source("Scripts/_functions/gene_annotation/add_gene_info_cols.R")

#'
ods_ss <- readRDS(snakemake@input$ods_ss)
dim(ods_ss)

ods_ns <- readRDS(snakemake@input$ods_ns)
dim(ods_ns)

par(mfrow = c(1,2))

#' ## Visualize
#' ### Size Factors
barplot(sort(sizeFactors(ods_ss)), main = "Strand Specific", xaxt = 'n', xlab = 'rank', ylab = 'Size Factors')
barplot(sort(sizeFactors(ods_ns)), main = "Non Strand Specific", xaxt = 'n', xlab = 'rank', ylab = 'Size Factors')

#' ### Aberrant samples
plotAberrantPerSample(ods_ss, main = 'Strand Specific')
plotAberrantPerSample(ods_ns, main = 'Non Strand Specific')


#' ### Batch correction
#+ heatmap, fig.height=8, fig.width=8
plotCountCorHeatmap(ods_ss, normalized=FALSE, rowCoFactor = "batch", main = 'Strand Specific before')
plotCountCorHeatmap(ods_ss, normalized=TRUE, rowCoFactor = "batch", main = 'Strand Specific after')

plotCountCorHeatmap(ods_ns, normalized=FALSE, rowCoFactor = "batch", main = 'Non Strand Specific before')
plotCountCorHeatmap(ods_ns, normalized=TRUE, rowCoFactor = "batch", main = 'Non Strand Specific after')


#' ## Results
# Strand Specific
res_ss <- results(ods_ss)
dim(res_ss)
res_ss[, STRANDED := 'Specific']

# Non Strand Specific
res_ns <- results(ods_ns)
dim(res_ns)
res_ns[, STRANDED := 'Non Specific']

# Combine them and add genetic information
res <- rbind(res_ss, res_ns)
res[, FC := 2^l2fc]

res <- add_all_gene_info(res, gene_name = 'geneID', dis_genes = F)

# Add sample annotation
sa <- fread(snakemake@config$SAMPLE_ANNOTATION)
res[, geneID := toupper(geneID)]
res <- left_join(res, sa[, .(RNA_ID, FIBROBLAST_ID, EXOME_ID, PEDIGREE, KNOWN_MUTATION,
                             CANDIDATE_GENE, BATCH)],
                 by = c("sampleID" = "RNA_ID")) %>% as.data.table

res[, TP := as.character(geneID == KNOWN_MUTATION)]
res[is.na(KNOWN_MUTATION), TP := "Unsolved"]

res[, tp_sample := as.character(any(geneID == KNOWN_MUTATION)), by = sampleID]
res[is.na(KNOWN_MUTATION), tp_sample := "Unsolved"]

#' ### How many samples with at least one gene
res[, uniqueN(sampleID), by = STRANDED]

#' ### Aberrant samples
ab_ss_samples <- res[STRANDED == 'Specific' & AberrantBySample > nrow(ods_ss)/1000, .N, by = .(sampleID, STRANDED)] %>% unique
ab_ns_samples <- res[STRANDED == 'Non Specific' & AberrantBySample > nrow(ods_ns)/1000, .N, by = .(sampleID,STRANDED)] %>% unique
ab_table <- rbind(ab_ss_samples, ab_ns_samples)
setorder(ab_table, N)
DT::datatable(left_join(ab_table, 
                        sa[,.(RNA_ID,BATCH,FIBROBLAST_ID,EXOME_ID,GENOME_ID,PROTEOME_ID,GENDER,DISEASE,KNOWN_MUTATION,PEDIGREE,RNA_PERSON,CANDIDATE_GENE)], 
                        by = c("sampleID" = "RNA_ID")), caption = "Aberrant samples", style = 'bootstrap')

#' ### Download results table
write.table(res, snakemake@output[['results']], quote = F, row.names = F, sep = "\t")
write.table(res, "/s/public_webshare/project/genetic_diagnosis/results/OUTRIDER_results.tsv", sep = "\t", quote = F, row.names = F)

#' [Download OTURIDER results table](https://i12g-gagneurweb.informatik.tu-muenchen.de/project/genetic_diagnosis/results/OUTRIDER_results.tsv)
DT::datatable(res, caption = "OUTRIDER results", style = 'bootstrap')


#' ### Visualize results
library(ggbeeswarm)
library(ggthemes)

#' Distributin of fold changes (FC)
#+ fig.width=10
ggplot(res, aes(FC)) + geom_histogram(aes(fill = TP)) + facet_wrap(~ STRANDED) + theme_bw() + scale_fill_fivethirtyeight()

ggplot(res[FC < 1 & !sampleID %in% ab_table$sampleID], aes(TP, FC)) + geom_boxplot(outlier.shape = NA) + 
    geom_beeswarm(aes(col = TP)) + facet_wrap(~ STRANDED) + theme_bw() + scale_color_fivethirtyeight()
ggplot(res[FC < 1 & !sampleID %in% ab_table$sampleID], aes(TP, -log10(padjust))) + geom_boxplot(outlier.shape = NA) + 
    geom_beeswarm(aes(col = TP)) + facet_wrap(~ STRANDED) + theme_bw() + scale_color_fivethirtyeight()

abt <- res[, .N, by = .(sampleID, tp_sample, STRANDED)]
setorder(abt, N)
ggplot(abt, aes(tp_sample, N)) + geom_boxplot(outlier.shape = NA) + geom_beeswarm(aes(col = tp_sample)) + 
    facet_wrap(~ STRANDED, scales = 'free') + labs(y = 'Number of outlier genes') + theme_bw(base_size = 14) + scale_color_fivethirtyeight()

#' Check the p adjusted of the true positives with more than one outlier
ggplot(res[tp_sample == T & AberrantBySample > 1], aes(sampleID, -log10(padjust))) + 
    geom_beeswarm(aes(col = TP)) + theme_bw(base_size=14) + scale_color_fivethirtyeight() + coord_flip()
