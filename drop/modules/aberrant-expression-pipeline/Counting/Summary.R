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

has_external <- !(all(ods@colData$GENE_COUNTS_FILE == "") || is.null(ods@colData$GENE_COUNTS_FILE))
if(has_external){
    ods@colData$isExternal <- as.factor(ods@colData$GENE_COUNTS_FILE != "")
}else{
    ods@colData$isExternal <- as.factor(FALSE)
}

# save ods with isExternal column
saveRDS(ods,snakemake@input$ods)

cnts_mtx_local <- counts(ods, normalized = F)[,!ods@colData$isExternal]
cnts_mtx <- counts(ods, normalized = F)

#' ## Number of samples:  
#' Local: `r sum(!as.logical(ods@colData$isExternal))`  
#' External: `r sum(as.logical(ods@colData$isExternal))`  
#' 
#' # Count Quality Control
#' 
#' Compare number of records vs. read counts  
#' `The Obtained Read Count Ratio` plot does not include external counts
#' because there are no raw reads to be counted.
#' 
bam_coverage <- fread(snakemake@input$bam_cov)
bam_coverage[, sampleID := as.character(sampleID)]
coverage_dt <- merge(bam_coverage,
                     data.table(sampleID = colnames(ods),
                                read_count = colSums(cnts_mtx),
                                isExternal = ods@colData$isExternal),
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

p_depth <- ggplot(coverage_dt, aes(x = count_rank, y = read_count, col= isExternal )) +
  geom_point(size = 3,show.legend = has_external) +
  theme_cowplot() +
  background_grid() +
  labs(title = "Obtained Read Counts", x="Sample Rank", y = "Reads Counted") +
  ylim(c(0,NA)) +
  scale_color_brewer(palette="Dark2")

p_frac <- ggplot(coverage_dt, aes(x = frac_rank, y = counted_frac, col = isExternal)) +
  geom_point(size = 3,show.legend = has_external) +
  theme_cowplot() +
  background_grid() +
  labs(title = "Obtained Read Count Ratio", x = "Sample Rank", 
       y = "Percent Reads Counted") +
  ylim(c(0,NA)) +
  scale_color_brewer(palette="Dark2")

#+ QC, fig.height=6, fig.width=12
plot_grid(p_depth, p_frac) 

p_sf <- ggplot(coverage_dt, aes(sf_rank, size_factors, col = isExternal)) +
  geom_point(size = 3,show.legend = has_external) +
  ylim(c(0,NA)) +
  theme_cowplot() +
  background_grid() +
  labs(title = 'Size Factors', x = 'Sample Rank', y = 'Size Factors') +
  scale_color_brewer(palette="Dark2")

p_sf_cov <- ggplot(coverage_dt, aes(read_count, size_factors, col = isExternal)) +
  geom_point(size = 3,show.legend = has_external) +
  ylim(c(0,NA)) +
  theme_cowplot() +
  background_grid() +
  labs(title = 'Size Factors vs. Read Counts',
       x = 'Read Counts', y = 'Size Factors') +
  scale_color_brewer(palette="Dark2")

#+ sizeFactors, fig.height=6, fig.width=12
plot_grid(p_sf, p_sf_cov)

#' # Filtering
#' **local**: A pre-filtered summary of counts using only the local (from BAM) counts. Omitted if no external counts  
#' **all**: A pre-filtered summary of counts using only the merged local (from BAM) and external counts  
#' **passed_FPKM**: Passes the user defined FPKM cutoff in at least 5% of genes  
#' **min 1 read**: minimum of 1 read expressed in 5% of genes  
#' **min 10 reads**: minimum of 10 reads expressed in 5% of genes  

quant <- .95

if(has_external){
    filter_mtx <- list(
      local = cnts_mtx_local,
      all = cnts_mtx,
      `passed FPKM` = cnts_mtx[rowData(ods)$passedFilter,],
      `min 1 read` = cnts_mtx[rowQuantiles(cnts_mtx, probs = quant) > 1, ],
      `min 10 reads` = cnts_mtx[rowQuantiles(cnts_mtx, probs = quant) > 10, ]
    )
    filter_dt <- lapply(names(filter_mtx), function(filter_name) {
      mtx <- filter_mtx[[filter_name]]
      data.table(gene_ID = rownames(mtx), median_counts = rowMeans(mtx), filter = filter_name)
    }) %>% rbindlist
    filter_dt[, filter := factor(filter, levels = c('local', 'all', 'passed FPKM', 'min 1 read', 'min 10 reads'))]
} else{
    filter_mtx <- list(
      all = cnts_mtx,
      `passed FPKM` = cnts_mtx[rowData(ods)$passedFilter,],
      `min 1 read` = cnts_mtx[rowQuantiles(cnts_mtx, probs = quant) > 1, ],
      `min 10 reads` = cnts_mtx[rowQuantiles(cnts_mtx, probs = quant) > 10, ]
    )
    filter_dt <- lapply(names(filter_mtx), function(filter_name) {
      mtx <- filter_mtx[[filter_name]]
      data.table(gene_ID = rownames(mtx), median_counts = rowMeans(mtx), filter = filter_name)
    }) %>% rbindlist
    filter_dt[, filter := factor(filter, levels = c('all', 'passed FPKM', 'min 1 read', 'min 10 reads'))]
}

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

#' ### Expressed Genes
exp_genes_cols <- c(Rank = "expressedGenesRank",`Expressed\ngenes` = "expressedGenes", 
                    `Union of\nexpressed genes` = "unionExpressedGenes", 
                    `Intersection of\nexpressed genes` = "intersectionExpressedGenes", 
                    `Genes passed\nfiltering` = "passedFilterGenes", `Is External` = "isExternal")

expressed_genes <- as.data.table(colData(ods)[,exp_genes_cols])
colnames(expressed_genes) <- names(exp_genes_cols)

#+ expressedGenes, fig.height=6, fig.width=8
plotExpressedGenes(ods) + 
  theme_cowplot() +
  background_grid(major = "y") +
  geom_point(data =melt(expressed_genes,id.vars = c("Rank","Is External")),
             aes(x = Rank, y = value, col = variable, shape = `Is External`),show.legend = has_external)

if(has_external){
    DT::datatable(expressed_genes[order(Rank)],rownames = F)
} else{
    DT::datatable(expressed_genes[order(Rank),-"Is External"],rownames = F)
}
