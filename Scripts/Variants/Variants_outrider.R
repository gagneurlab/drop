#'---
#' title: Variants associated to OUTRIDER
#' author: vyepez
#' wb:
#'  input:
#'   - done_file: '`sm config["PROC_RESULTS"] + "/process_vcf/vars_outrider/outrider_all_vars.Rds"`'
#'  output:
#'  threads: 10
#' output:
#'   html_document:
#'    code_folding: show
#'    code_download: TRUE
#'---

#+ echo=F
saveRDS(snakemake, "tmp/outrider_vars_plots.snakemake")
# snakemake <- readRDS("tmp/outrider_vars_plots.snakemake")

suppressPackageStartupMessages({
  library(ggplot2)
  library(data.table)
  library(dplyr)
  library(cowplot)
  library(ggpval)
  library(BiocParallel)
})


var_type = 'stop'

PT <- readRDS(snakemake@input$done_file)

for(var_type in c("missense", "synonymous", "splice", "frame-shift", "stop")){
  print(var_type)
  # pt <- readRDS(paste0("/s/project/genetic_diagnosis/processed_results/process_vcf/vars_outrider/", var_type, ".Rds"))
  
  # has_col <- paste0("has_", var_type)
  # has_rare_col <- paste0("has_rare_", var_type)
  
  # print(paste0("Number of rare ", var_type, " variants"))
  # table(pt[, get(has_rare_col)]) %>% print
  
  # print(paste0("Number of rare ", var_type, " variants with outliers"))
  # table(pt[get(has_rare_col), aberrant]) %>% print
  
  # table(pt[has_rare_missense == T, padjust < .1])
  # table(pt[has_rare_missense == T, padjust < .3])
  # table(pt[has_rare_missense == T, padjust < .5])
  # table(pt[has_rare_missense == T, padjust < 1])
  # View(pt[has_rare_missense == T])
  
  # g <- ggplot(pt, aes(get(has_col), l2fc)) + geom_boxplot() + theme_bw(base_size = 14) + labs(x = has_col)
  # plot(g)
  # g2 <- ggplot(pt, aes(get(has_rare_col), l2fc)) + geom_boxplot() + theme_bw(base_size = 14) + labs(x = has_rare_col)
  # plot(g2)
  
  # g3 <- ggplot(pt[get(has_col) == T], aes(get(has_rare_col), CADD_raw)) + geom_boxplot(aes(fill = aberrant)) + theme_bw(base_size = 14) + labs(x = has_rare_col)
  # g3
}

PT[, padj_rank := NULL]
PT = unique(PT)
pt_out <- PT[aberrant == T]
dim(pt_out)
mt <- melt(pt_out[,.(geneID, sampleID, has_missense, has_rare_missense, has_stop, has_rare_stop, `has_frame-shift`, `has_rare_frame-shift`, 
               has_synonymous, has_rare_synonymous, has_splice, has_rare_splice)], id.vars = c("geneID", "sampleID"))
mt[, variable := gsub("has_", "", variable)]
mt = mt[value == T]
mt_rare <- mt[grep('rare', variable)]
# ggplot(mt_rare, aes(variable)) + geom_bar(aes(y = ..count..)) + labs(y = 'Number of outliers', x = 'Variant type') + coord_flip()
# ggplot(mt[!grep('rare', variable)], aes(variable)) + geom_bar(aes(y = ..count..)) + labs(y = 'Number of outliers', x = 'Variant type') + coord_flip()
# pt_out[, c('rawcounts', 'normcounts', 'meanCorrected',  'theta', 'aberrant', 'AberrantBySample', 'AberrantByGene') := NULL]
# pt_out[has_rare_stop == T & KNOWN_MUTATION != geneID]
# pt_out[has_rare_stop == T & KNOWN_MUTATION == geneID]
# pt_out[has_rare_stop == T & is.na(KNOWN_MUTATION)]
