#'---
#' title: "Counts Summary: `r paste(snakemake@wildcards$dataset, snakemake@wildcards$annotation, sep = '--')`"
#' author: 
#' wb:
#'  log:
#'   - snakemake: '`sm str(tmp_dir / "AE" / "{annotation}" / "{dataset}" / "count_summary.Rds")`'
#'  input: 
#'    - ods: '`sm cfg.getProcessedResultsDir() +
#'            "/aberrant_expression/{annotation}/outrider/{dataset}/ods_unfitted.Rds"`'
#'    - bam_cov: '`sm rules.aberrantExpression_mergeBamStats.output`'
#'  output:
#'   - wBhtml: '`sm config["htmlOutputPath"] +
#'              "/AberrantExpression/Counting/{annotation}/Summary_{dataset}.html"`'
#'  type: noindex
#' output:
#'  html_document:
#'   code_folding: hide
#'   code_download: TRUE
#'---

saveRDS(snakemake, snakemake@log$snakemake)

suppressPackageStartupMessages({
  library(OUTRIDER)
  library(SummarizedExperiment)
  library(GenomicAlignments)
  library(ggplot2)
  library(ggthemes)
  library(cowplot)
  library(data.table)
  library(tidyr)
})

ods <- readRDS(snakemake@input$ods)
cnts_mtx <- counts(ods, normalized = F)

#' Number of samples: `r ncol(ods)`
#' 
#' # Count Quality Control
#' 
#' Compare number of records vs. read counts
#' 
bam_coverage <- fread(snakemake@input$bam_cov)
bam_coverage[, sampleID := as.character(sampleID)]
coverage_dt <- merge(bam_coverage,
                     data.table(sampleID = colnames(ods),
                                read_count = colSums(cnts_mtx)),
                     by = "sampleID", sort = FALSE)
# read count
setorder(coverage_dt, read_count)
coverage_dt[, count_rank := .I]
# ratio
coverage_dt[, counted_frac := read_count/record_count]
setorder(coverage_dt, counted_frac)
coverage_dt[, frac_rank := .I]

# size factors 
ods <- estimateSizeFactors(ods)
coverage_dt[, size_factors := sizeFactors(ods)]
setorder(coverage_dt, size_factors)
coverage_dt[, sf_rank := 1:.N]

p_depth <- ggplot(coverage_dt, aes(count_rank, read_count)) +
  geom_point() +
  theme_cowplot() +
  background_grid() +
  labs(title = "Obtained Read Counts", x="Sample Rank", y = "Reads Counted") +
  ylim(c(0,NA))

p_frac <- ggplot(coverage_dt, aes(frac_rank, counted_frac)) +
  geom_point() +
  theme_cowplot() +
  background_grid() +
  labs(title = "Obtained Read Count Ratio", x = "Sample Rank", 
       y = "Percent Reads Counted") +
  ylim(c(0,NA))

#+ QC, fig.height=6, fig.width=12
plot_grid(p_depth, p_frac)

p_sf <- ggplot(coverage_dt, aes(sf_rank, size_factors)) +
  geom_point() +
  ylim(c(0,NA)) +
  theme_cowplot() +
  background_grid() +
  labs(title = 'Size Factors', x = 'Sample Rank', y = 'Size Factors')

p_sf_cov <- ggplot(coverage_dt, aes(read_count, size_factors)) +
  geom_point() +
  ylim(c(0,NA)) +
  theme_cowplot() +
  background_grid() +
  labs(title = 'Size Factors vs. Read Counts',
       x = 'Read Counts', y = 'Size Factors')

#+ sizeFactors, fig.height=6, fig.width=12
plot_grid(p_sf, p_sf_cov)

#' # Filtering
quant <- .95
filter_mtx <- list(
  all = cnts_mtx,
  passed_FPKM = cnts_mtx[rowData(ods)$passedFilter,],
  min_1 = cnts_mtx[rowQuantiles(cnts_mtx, probs = quant) > 1, ],
  min_10 = cnts_mtx[rowQuantiles(cnts_mtx, probs = quant) > 10, ]
)
filter_dt <- lapply(names(filter_mtx), function(filter_name) {
  mtx <- filter_mtx[[filter_name]]
  data.table(gene_ID = rownames(mtx), median_counts = rowMeans(mtx), filter = filter_name)
}) %>% rbindlist
filter_dt[, filter := factor(filter, levels = c('all', 'passed_FPKM', 'min_1', 'min_10'))]

binwidth <- .2
p_hist <- ggplot(filter_dt, aes(x = median_counts, fill = filter)) +
  geom_histogram(binwidth = binwidth) +
  scale_x_log10() +
  facet_wrap(.~filter) +
  labs(x = "Mean counts per gene", y = "Frequency", title = 'Mean Count Distribution') +
  guides(col = guide_legend(title = NULL)) +
  scale_fill_brewer(palette = "Paired") +
  theme_cowplot() +
  theme(legend.position = "none")

p_dens <- ggplot(filter_dt, aes(x = median_counts, col = filter)) +
  geom_density(aes(y=binwidth * ..count..), size = 1.2) +
  scale_x_log10() +
  labs(x = "Mean counts per gene", y = "Frequency") +
  guides(col = guide_legend(title = NULL)) +
  scale_color_brewer(palette = "Paired") +
  theme_cowplot() +
  theme(legend.position = "top",
        legend.justification="center",
        legend.background = element_rect(color = NA))

#+ meanCounts, fig.height=6, fig.width=12
plot_grid(p_hist, p_dens)

#+ expressedGenes, fig.height=6, fig.width=8
plotExpressedGenes(ods) +
  theme_cowplot() +
  background_grid(major = "y")

expressed_genes <- as.data.table(colData(ods))
expressed_genes <- expressed_genes[, .(expressedGenes, unionExpressedGenes,
                                       intersectionExpressedGenes, passedFilterGenes,
                                       expressedGenesRank)]

#+echo=F
rank_1 <- expressed_genes[expressedGenesRank == 1]
#' **Rank 1:**
#' `r as.character(rank_1$expressedGenes)` expressed genes
#+echo=F
rank_n <- expressed_genes[expressedGenesRank == .N]
#' **Rank `r rank_n$expressedGenesRank`:**  
#' `r as.character(rank_n$expressedGenes)` expressed genes  
#' `r as.character(rank_n$unionExpressedGenes)` expressed genes (union)  
#' `r as.character(rank_n$intersectionExpressedGenes)` expressed genes (intersection)  
#' `r as.character(rank_n$passedFilterGenes)` genes passed the filter
