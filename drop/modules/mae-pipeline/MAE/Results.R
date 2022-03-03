#'---
#' title: MAE Results table
#' author: vyepez
#' wb:
#'  log:
#'   - snakemake: '`sm str(tmp_dir / "MAE" / "{dataset}" / "{annotation}_results.Rds")`'
#'  params:
#'   - allelicRatioCutoff: '`sm cfg.MAE.get("allelicRatioCutoff")`'
#'   - padjCutoff: '`sm cfg.MAE.get("padjCutoff")`'
#'   - maxCohortFreq: '`sm cfg.MAE.get("maxVarFreqCohort")`'
#'  input:
#'   - mae_res: '`sm lambda w: expand(cfg.getProcessedResultsDir() + 
#'                "/mae/samples/{id}_res.Rds", id=cfg.MAE.getMaeByGroup({w.dataset}))`'
#'   - gene_name_mapping: '`sm cfg.getProcessedDataDir() +
#'                          "/mae/gene_name_mapping_{annotation}.tsv"`'
#'   - input_sample_params: '`sm cfg.getProcessedDataDir() + "/mae/params/results/{dataset}_resultParams.csv" `'
#'  output:
#'   - res_all: '`sm cfg.getProcessedResultsDir() + 
#'                "/mae/{dataset}/MAE_results_all_{annotation}.tsv.gz"`' 
#'   - res_signif: '`sm cfg.getProcessedResultsDir() + 
#'                   "/mae/{dataset}/MAE_results_{annotation}.tsv"`'
#'   - res_signif_rare: '`sm cfg.getProcessedResultsDir() + 
#'                   "/mae/{dataset}/MAE_results_{annotation}_rare.tsv"`'
#'   - wBhtml: '`sm config["htmlOutputPath"] +
#'               "/MonoallelicExpression/{dataset}--{annotation}_results.html"`'
#'  type: noindex
#'---

#+ echo=F
saveRDS(snakemake, snakemake@log$snakemake)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(tidyr)
  library(GenomicRanges)
  library(SummarizedExperiment)
  library(R.utils)
})

# Read all MAE results files
rmae <- lapply(snakemake@input$mae_res, function(m){
  rt <- readRDS(m)
  # force consistant UCSC chromosome style
  rt <- rt[!grepl("chr",contig),contig:= paste0("chr",contig)]
  return(rt)
}) %>% rbindlist()

# re-factor contig
rmae$contig <- factor(rmae$contig)

# Convert results into GRanges
rmae_ranges <- GRanges(seqnames = rmae$contig, 
                       IRanges(start = rmae$position, end = rmae$position), strand = '*')

# Read annotation and convert into GRanges
gene_annot_dt <- fread(snakemake@input$gene_name_mapping)
gene_annot_ranges <- GRanges(seqnames = gene_annot_dt$seqnames, 
                             IRanges(start = gene_annot_dt$start, end = gene_annot_dt$end), 
                             strand = gene_annot_dt$strand)
gene_annot_ranges <- keepStandardChromosomes(gene_annot_ranges, pruning.mode = 'coarse')

# Keep the chr style of the annotation in case the results contain different styles
seqlevelsStyle(rmae_ranges) <- seqlevelsStyle(gene_annot_ranges)

# Overlap results and annotation
fo <- findOverlaps(rmae_ranges, gene_annot_ranges)

# Add the gene names
res_annot <- cbind(rmae[from(fo), ],  gene_annot_dt[to(fo), .(gene_name, gene_type)])

# Prioritize protein coding genes
res_annot <- rbind(res_annot[gene_type == 'protein_coding'], 
                   res_annot[gene_type != 'protein_coding'])

# Write all the other genes in another column
res_annot[, aux := paste(contig, position, sep = "-")]
rvar <- unique(res_annot[, .(aux, gene_name)])
rvar[, N := 1:.N, by = aux]

r_other <- rvar[N > 1, .(other_names = paste(gene_name, collapse = ',')), by = aux]
res <- merge(res_annot, r_other, by = 'aux', sort = FALSE, all.x = TRUE) 
res[, c('aux') := NULL]
res <- res[, .SD[1], by = .(ID, contig, position)]

# Bring gene_name column front
res <- cbind(res[, .(gene_name)], res[, -"gene_name"])

# Calculate variant frequency within cohort 
maxCohortFreq <- snakemake@params$maxCohortFreq
res[, N_var := .N, by = .(gene_name, contig, position)]
res[, cohort_freq := round(N_var / uniqueN(ID), 3)]

res[, rare := (rare | is.na(rare)) & cohort_freq <= maxCohortFreq] 

# Add significance columns
allelicRatioCutoff <- snakemake@params$allelicRatioCutoff
res[, MAE := padj <= snakemake@params$padjCutoff &
      (altRatio >= allelicRatioCutoff | altRatio <= (1-allelicRatioCutoff)) 
    ] 
res[, MAE_ALT := MAE == TRUE & altRatio >= allelicRatioCutoff]
#'
#' Number of samples: `r uniqueN(res$ID)`
#'
#' Number of genes: `r uniqueN(res$gene_name)`
#'
#' Number of samples with significant MAE for alternative events: `r uniqueN(res[MAE_ALT == TRUE, ID])`

#+echo=F

# Save full results zipped
res[, altRatio := round(altRatio, 3)]
fwrite(res, snakemake@output$res_all, sep = '\t', 
       row.names = F, quote = F, compress = 'gzip')

# Save significant results
fwrite(res[MAE_ALT == TRUE], snakemake@output$res_signif, 
       sep = '\t', row.names = F, quote = F)

# Save significant results
fwrite(res[MAE_ALT == TRUE & rare == TRUE], snakemake@output$res_signif_rare, 
       sep = '\t', row.names = F, quote = F)


# Add columns for plot
res[, N := .N, by = ID]
res[MAE == TRUE, N_MAE := .N, by = ID]
res[MAE == TRUE & MAE_ALT == FALSE, N_MAE_REF := .N, by = ID]
res[MAE_ALT == TRUE, N_MAE_ALT := .N, by = ID]
res[MAE == TRUE & MAE_ALT == FALSE & rare == TRUE, N_MAE_REF_RARE := .N, by = ID]
res[MAE_ALT == TRUE & rare == TRUE, N_MAE_ALT_RARE := .N, by = ID]

rd <- unique(res[,.(ID, N, N_MAE, N_MAE_REF, N_MAE_ALT, N_MAE_REF_RARE, N_MAE_ALT_RARE)])
melt_dt <- melt(rd, id.vars = 'ID')
melt_dt[variable == 'N', variable := '>10 counts']
melt_dt[variable == 'N_MAE', variable := '+MAE']
melt_dt[variable == 'N_MAE_REF', variable := '+MAE for\nREF']
melt_dt[variable == 'N_MAE_ALT', variable := '+MAE for\nALT']
melt_dt[variable == 'N_MAE_REF_RARE', variable := '+MAE for REF\n& rare']
melt_dt[variable == 'N_MAE_ALT_RARE', variable := '+MAE for ALT\n& rare']

#' 
#' ## Cascade plot 
ggplot(melt_dt, aes(variable, value)) + geom_boxplot() +
  scale_y_log10() + theme_bw(base_size = 14) +
  labs(y = 'Heterozygous SNVs per patient', x = '') +
    annotation_logticks(sides = "l")

#'
#' ## Variant Frequency within Cohort Histogram
ggplot(unique(res[,cohort_freq,by =.(gene_name, contig, position)]),aes(x = cohort_freq)) + geom_histogram( binwidth = 0.02)  +
  geom_vline(xintercept = maxCohortFreq, col = "red") +
  xlab("Variant frequency in cohort") + ylab("Count")

#' Median of each category
DT::datatable(melt_dt[, .(median = median(value, na.rm = T)), by = variable])

#' 
#' ## MAE Results table
DT::datatable(
  head(res[MAE_ALT == TRUE], 1000),
  caption = 'MAE results (up to 1,000 rows shown)',
  options=list(scrollX=TRUE),
  filter = 'top'
)

