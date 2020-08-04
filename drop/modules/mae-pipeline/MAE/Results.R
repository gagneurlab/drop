#'---
#' title: MAE Results table
#' author: vyepez
#' wb:
#'  log:
#'   - snakemake: '`sm str(tmp_dir / "MAE" / "{dataset}" / "{annotation}_results.Rds")`'
#'  params:
#'   - allelicRatioCutoff: '`sm cfg.MAE.get("allelicRatioCutoff")`'
#'   - padjCutoff: '`sm cfg.MAE.get("padjCutoff")`'
#'  input:
#'   - mae_res: '`sm lambda w: expand(cfg.getProcessedResultsDir() + 
#'                "/mae/samples/{id}_res.Rds", id=cfg.MAE.getMaeByGroup({w.dataset}))`'
#'   - gene_name_mapping: '`sm cfg.getProcessedDataDir() +
#'                          "/mae/gene_name_mapping_{annotation}.tsv"`'
#'  output:
#'   - res_all: '`sm cfg.getProcessedResultsDir() + 
#'                "/mae/{dataset}/MAE_results_all_{annotation}.tsv.gz"`' 
#'   - res_signif: '`sm cfg.getProcessedResultsDir() + 
#'                   "/mae/{dataset}/MAE_results_{annotation}.tsv"`'
#'   - wBhtml: '`sm config["htmlOutputPath"] +
#'               "/MAE/{dataset}--{annotation}_results.html"`'
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
  return(rt)
}) %>% rbindlist()


# Add gene names
gene_annot_dt <- fread(snakemake@input$gene_name_mapping)
gene_annot_dt <- gene_annot_dt[seqnames %in% paste0('chr', c(1:22, 'X'))]

# Subtract the genomic ranges from the annotation and results and overlap them
gene_annot_ranges <- GRanges(seqnames = gene_annot_dt$seqnames, 
                             IRanges(start = gene_annot_dt$start, end = gene_annot_dt$end), 
                             strand = gene_annot_dt$strand)
rmae_ranges <- GRanges(seqnames = rmae$contig, 
                       IRanges(start = rmae$position, end = rmae$position), strand = '*')

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

# Bring gene_name column front
res <- cbind(res[, .(gene_name)], res[, -"gene_name"])

#'
#' Number of samples: `r uniqueN(res$ID)`
#'
#' Number of genes: `r uniqueN(res$gene_name)`

# Subset for significant events
allelicRatioCutoff <- snakemake@params$allelicRatioCutoff
res[, MAE := padj <= snakemake@params$padjCutoff &
      (altRatio >= allelicRatioCutoff | altRatio <= (1-allelicRatioCutoff))] 
res[, MAE_ALT := MAE == TRUE & altRatio >= allelicRatioCutoff]

#' Number of samples with significant MAE for alternative events: `r uniqueN(res[MAE_ALT == TRUE, ID])`

### Save the results
# Save full results zipped
fwrite(res, snakemake@output$res_all, sep = '\t', 
       row.names = F, quote = F, compress = 'gzip')

# Save significant results
fwrite(res[MAE_ALT == TRUE], snakemake@output$res_signif, 
       sep = '\t', row.names = F, quote = F)


#+echo=F
res[, N := .N, by = ID]
res[MAE == TRUE, N_MAE := .N, by = ID]
res[MAE_ALT == TRUE, N_MAE_ALT := .N, by = ID]
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

#' Median of each category
DT::datatable(melt_dt[, .(median = median(value, na.rm = T)), by = variable])

#' 
#' ## Results table
DT::datatable(res[MAE_ALT == TRUE], filter = 'top')

