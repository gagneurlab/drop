#'---
#' title: OUTRIDER Summary
#' author: mumichae, vyepez
#' wb:
#'  input:
#'   - ods: '`sm config["PROC_RESULTS"] + "/{annotation}/outrider/{dataset}/ods.Rds"`'
#'   - functions: "Scripts/_functions/gene_annotation/add_gene_info_cols.R"
#'   - results: '`sm config["PROC_RESULTS"] + "/{annotation}/outrider/{dataset}/OUTRIDER_results.tsv"`'
#'   - results_public: "/s/public_webshare/project/genetic_diagnosis/results/{annotation}/OUTRIDER_results_{dataset}.tsv"
#'  output:
#'   - wBhtml: "Output/html/AberrantExpression/Outrider/{annotation}/OutriderSummary_{dataset}.html"
#'  type: noindex
#'---

#+ echo=F
saveRDS(snakemake, "tmp/outrider_summary.snakemake")
# snakemake <- readRDS("tmp/outrider_summary.snakemake")

suppressPackageStartupMessages({
    library(OUTRIDER)
    library(SummarizedExperiment)
    library(ggplot2)
    library(cowplot)
    library(data.table)
    library(dplyr)
    library(ggbeeswarm)
    library(ggthemes)
})

source(snakemake@input$functions)

#' ## Read ods object
ods <- readRDS(snakemake@input$ods)
# Number of samples and genes
dim(ods)

#' ## Visualize
#' ### Parameters
barplot(sort(sizeFactors(ods)), main = paste('Size Factors (', snakemake@wildcards$dataset, ')'), xaxt = 'n', xlab = 'rank', ylab = 'Size Factors')
plotEncDimSearch(ods)

#' ### Aberrant samples
plotAberrantPerSample(ods, main = snakemake@wildcards$dataset)


devtools::load_all("../OUTRIDER/")

#' ### Batch correction
#+ heatmap, fig.height=8, fig.width=8
plotCountCorHeatmap(ods, normalized = FALSE, rowGroups = NA, colGroups = c("GENDER", "BATCH", "TISSUE", "GROWTH_MEDIUM"), 
                    main = paste('Raw Counts (', snakemake@wildcards$dataset, ')'))
plotCountCorHeatmap(ods, normalized = TRUE, rowCoFactor = NA, colGroups = c("GENDER", "BATCH", "TISSUE", "GROWTH_MEDIUM"),
                    main = paste('Normalized Counts (', snakemake@wildcards$dataset, ')'))

#' ## Results
res <- fread(snakemake@input$results)
sa <- fread(snakemake@config$SAMPLE_ANNOTATION)

#' ### How many samples with at least one gene
res[, uniqueN(sampleID)]

#' ### Aberrant samples
if (nrow(res) > 0) {
    ab_table <- res[AberrantBySample > nrow(ods)/1000, .N, by = .(sampleID)] %>% unique
    setorder(ab_table, N)
    DT::datatable(left_join(ab_table, 
                            sa[,.(RNA_ID,BATCH,FIBROBLAST_ID,EXOME_ID,GENOME_ID,PROTEOME_ID,GENDER,DISEASE,KNOWN_MUTATION,PEDIGREE,RNA_PERSON,CANDIDATE_GENE,COMMENT)], 
                            by = c("sampleID" = "RNA_ID")), caption = "Aberrant samples", style = 'bootstrap')
} else {
    print('no results')
}

#' ### Download results table
results_link <- paste0('https://i12g-gagneurweb.informatik.tu-muenchen.de/project/genetic_diagnosis/results/', snakemake@wildcards$annotation,'/OUTRIDER_results_', snakemake@wildcards$dataset, '.tsv')
#' [Download OUTRIDER results table](`r results_link`)
DT::datatable(res, caption = "OUTRIDER results", style = 'bootstrap', filter = 'top')

#' ### Visualize Results
#' Distribution of fold changes (FC)
#+ fig.width=10
# ggplot(res, aes(FC)) + geom_histogram(aes(fill = TP)) + scale_fill_fivethirtyeight()

# if(exists("ab_table")) {
#     ggplot(res[FC < 1 & !sampleID %in% ab_table$sampleID], aes(TP, FC)) + geom_boxplot(outlier.shape = NA) + 
#        geom_beeswarm(aes(col = TP)) + scale_color_fivethirtyeight()
#    ggplot(res[FC < 1 & !sampleID %in% ab_table$sampleID], aes(TP, -log10(padjust))) + geom_boxplot(outlier.shape = NA) + 
#        geom_beeswarm(aes(col = TP)) + scale_color_fivethirtyeight() 
#}

# abt <- res[, .N, by = .(sampleID, tp_sample)]
# setorder(abt, N)
#ggplot(abt, aes(tp_sample, N)) + 
#    geom_boxplot(outlier.shape = NA) +
#    geom_beeswarm(aes(col = tp_sample)) + 
#    labs(y = 'Number of outlier genes') +
#    scale_color_fivethirtyeight()

#' Check the p adjusted of the true positives with more than one outlier
# ggplot(res[tp_sample == T & AberrantBySample > 1], aes(sampleID, -log10(padjust))) + 
#    geom_beeswarm(aes(col = TP)) +
#    scale_color_fivethirtyeight() + coord_flip()




