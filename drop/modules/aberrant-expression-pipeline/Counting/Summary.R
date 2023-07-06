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

has_external <- any(as.logical(colData(ods)$isExternal))
cnts_mtx_local <- counts(ods, normalized = F)[,!as.logical(ods@colData$isExternal),drop=FALSE]
cnts_mtx <- counts(ods, normalized = F)

#' ## Number of samples:  
#' Local: `r sum(!as.logical(ods@colData$isExternal))`  
#' External: `r sum(as.logical(ods@colData$isExternal))`  
#' 
#' # Count Quality Control
#' 
#' Compares total mapped vs counted reads.  
#' The `Mapped vs Counted Reads` plot does not include external counts.  
#' Consider removing samples with too low or too high size factors.
#'  
bam_coverage <- fread(snakemake@input$bam_cov)
bam_coverage[, RNA_ID := as.character(sampleID)]
bam_coverage[, sampleID := NULL]
setnames(bam_coverage, 'record_count', 'total_count')
coverage_dt <- merge(data.table(RNA_ID = colnames(ods),
                                read_count = colSums(cnts_mtx),
                                isExternal = ods@colData$isExternal),
                     bam_coverage,
                     by = "RNA_ID", sort = FALSE)
# read counts
coverage_dt[, count_rank := rank(read_count)]

# size factors 
ods <- estimateSizeFactors(ods)
coverage_dt[, size_factors := round(sizeFactors(ods), 3)]
coverage_dt[, sf_rank := rank(size_factors)]

# Show this plot only if there are external samples, as the next plot contains the same info
if(has_external == T){
  p_depth <- ggplot(coverage_dt, aes(x = count_rank, y = read_count, col = isExternal)) +
    geom_point(size = 3,show.legend = has_external) +
    theme_cowplot() + background_grid() +
    labs(title = "Obtained Read Counts", x="Sample Rank", y = "Counted Reads") +
    ylim(c(0,NA)) +
    scale_color_brewer(palette="Dark2")
  p_depth
}


p_comp <- ggplot(coverage_dt[isExternal == F], aes(x = total_count, y = read_count)) +
  geom_point(size = 3, show.legend = has_external, color = "#1B9E77") +
  theme_cowplot() + background_grid() +
  labs(title = "Total mapped vs. Counted Reads", x = "Mapped Reads", y = "Counted Reads") +
  xlim(c(0,NA)) + ylim(c(0,NA))
p_comp
# ggExtra::ggMarginal(p_comp, type = "histogram") # could be a nice add-on

p_sf <- ggplot(coverage_dt, aes(sf_rank, size_factors, col = isExternal)) +
  geom_point(size = 3,show.legend = has_external) +
  ylim(c(0,NA)) +
  theme_cowplot() + background_grid() +
  labs(title = 'Size Factors', x = 'Sample Rank', y = 'Size Factors') +
  scale_color_brewer(palette="Dark2")

p_sf

setnames(coverage_dt, old = c('total_count', 'read_count', 'size_factors'),
         new = c('Reads Mapped', 'Reads Counted', 'Size Factors'))
DT::datatable(coverage_dt[, .(RNA_ID, `Reads Mapped`, `Reads Counted`, `Size Factors`)][order(`Reads Mapped`)],
              caption = 'Reads summary statistics')

#' # Filtering
#' **local**: A pre-filtered summary of counts using only the local (from BAM) counts. Omitted if no external counts  
#' **all**: A pre-filtered summary of counts using only the merged local (from BAM) and external counts  
#' **passed FPKM**: Passes the user defined FPKM cutoff in at least 5% of genes  
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

expressed_genes <- as.data.table(colData(ods)[,exp_genes_cols], keep.rownames = TRUE)
colnames(expressed_genes) <- c('RNA_ID', names(exp_genes_cols))

#+ expressedGenes, fig.height=6, fig.width=8
plotExpressedGenes(ods) + 
  theme_cowplot() +
  background_grid(major = "y") +
  geom_point(data = melt(expressed_genes, id.vars = c("RNA_ID", "Rank", "Is External")),
             aes(Rank, value, col = variable, shape = `Is External`), 
             show.legend = has_external)

if(has_external){
    DT::datatable(expressed_genes[order(Rank)],rownames = F)
} else{
    DT::datatable(expressed_genes[order(Rank),-"Is External"],rownames = F)
}

#' **Considerations:**
#' The verifying of the samples sex is performed by comparing the expression levels of 
#' the genes XIST and UTY. In general, females should have high XIST and low UTY expression,
#' and viceversa for males. For it to work, the sample annotation must have a column called 'SEX',
#' with values male/female. If some other values exist, e.g., unknown, they will be matched. 
#' The prediction is performed via a linear discriminant analysis on the log2 counts.

# Get sex column and proceed if it exists
sex_idx <- which('SEX' == toupper(colnames(colData(ods))))
if(isEmpty(sex_idx)){
  print('Sex column not found in sample annotation')
} else{
  
  # Verify if both XIST and UTY were counted
  xist_id <- 'XIST'
  uty_id <- 'UTY'
  
  if(grepl('ENSG0', rownames(ods)[1])){
    xist_id <- 'ENSG00000229807'
    uty_id <- 'ENSG00000183878'
  }
  xist_idx <- grep(xist_id, rownames(ods))
  uty_idx <- grep(uty_id, rownames(ods))
  
  if(isEmpty(xist_idx) | isEmpty(uty_idx)){
    print('Either XIST or UTY is not expressed')
  }else{
    
    sex_counts <- counts(ods)[c(xist_idx, uty_idx), ]
    sex_dt <- data.table(sampleID = colnames(ods), 
                         XIST = counts(ods)[xist_idx,], 
                         UTY = counts(ods)[uty_idx,])
    sex_dt <- merge(sex_dt, as.data.table(colData(ods))[,c(1, sex_idx), with = F], sort = F)
    colnames(sex_dt) <- toupper(colnames(sex_dt))
    sex_dt[, SEX := tolower(SEX)]
    sex_dt[is.na(SEX), SEX := '']
    
    # Train only in male/female in case there are other values
    train_dt <- sex_dt[SEX %in% c('f', 'm', 'female', 'male')]
    
    library("MASS")
    lda <- lda(SEX ~ log2(XIST+1) + log2(UTY+1), data = train_dt)
    
    sex_dt[, predicted_sex := predict(lda, sex_dt)$class]
    sex_dt[, match_sex := SEX == predicted_sex]
    table(sex_dt[, .(SEX, predicted_sex)])
    
    g <- ggplot(sex_dt, aes(XIST+1, UTY+1)) + 
      geom_point(aes(col = SEX, shape = predicted_sex, alpha = match_sex)) + 
      scale_x_log10(limits = c(1,NA)) + scale_y_log10(limits = c(1,NA)) +
      scale_alpha_manual(values = c(1, .1)) + 
      theme_cowplot() + background_grid(major = 'xy', minor = 'xy') + 
      annotation_logticks(sides = 'bl') + 
      labs(color = 'Sex', shape = 'Predicted sex', alpha = 'Matches sex')
    plot(g)
    
    DT::datatable(sex_dt[match_sex == F], caption = 'Sex mismatches')
    
    # Write table
    fwrite(sex_dt, gsub('ods_unfitted.Rds', 'xist_uty.tsv', snakemake@input$ods), 
           sep = '\t', quote = F)
  }
}
