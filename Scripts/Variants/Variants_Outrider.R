#'---
#' title: Variants associated to OUTRIDER
#' author: vyepez
#' wb:
#'  input:
#'   - var_dt: '`sm config["PROC_RESULTS"] + "/process_vcf/variant_dt.Rds"`'
#'   - res_fib: '`sm config["PROC_RESULTS"] + "/v29_overlap/outrider/fib_all/OUTRIDER_results_all.Rds"`'
#'  output:
#'  threads: 10
#' output:
#'   html_document:
#'    code_folding: show
#'    code_download: TRUE
#'---

#+ echo=F
saveRDS(snakemake, "tmp/outrider_vars.snakemake")
# snakemake <- readRDS("tmp/outrider_vars.snakemake")

suppressPackageStartupMessages({
  library(ggplot2)
  library(data.table)
  library(dplyr)
  library(cowplot)
  library(OUTRIDER)
  library(ggpval)
  library(BiocParallel)
})


#' Read results and remove samples without exomes
sa <- fread(snakemake@config$SAMPLE_ANNOTATION)
res_fib <- readRDS(snakemake@input$res_fib)   # TODO: change to fib, not fib_all
res_fib <- readRDS("/s/project/genetic_diagnosis/processed_results/v29_overlap/outrider/fib_all/OUTRIDER_results_all.Rds")
res_fib <- left_join(res_fib, sa[GROWTH_MEDIUM == 'GLU' & is.na(TRANSDUCED_GENE), .SD[1], by = 'EXOME_ID'][,.(RNA_ID, EXOME_ID, KNOWN_MUTATION)], by = c("sampleID" = "RNA_ID")) %>% as.data.table
res_fib[, geneID := toupper(geneID)]
res_fib = res_fib[!is.na(EXOME_ID)]

register(MulticoreParam(snakemake@threads))

#'+results=F
bplapply(c("missense", "synonymous", "splice", "unstop", "frame-shift", "unstart", "stop", "stop_retain", "coding"), 
    function(var_type){
       mt <- readRDS(paste0("/s/project/genetic_diagnosis/processed_results/process_vcf/unique_vars_gene/", var_type, "_unique_variant_dt.Rds"))
       exomes_common <- intersect(res_fib$EXOME_ID, mt$EXOME_ID)
       length(exomes_common)
       # Join only samples that are both on OUTRIDER and with variants
       pt <- left_join(res_fib[EXOME_ID %in% exomes_common], mt, by = c("EXOME_ID", "geneID")) %>% as.data.table
       
       has_col <- paste0("has_", var_type)
       has_rare_col <- paste0("has_rare_", var_type)
       set(pt, j = paste0("has_", var_type), value = !is.na(pt$pos))
       set(pt, j = has_rare_col, value = !is.na(pt$pos) & pt$MAX_AF < 0.001)
       
       print(paste0("Number of ", var_type, " variants"))
       table(pt[, get(has_rare_col)]) %>% print
       
       print(paste0("Number of ", var_type, " variants with outliers"))
       table(pt[get(has_rare_col), aberrant]) %>% print
       
       # table(pt[has_rare_missense == T, padjust < .1])
       # table(pt[has_rare_missense == T, padjust < .3])
       # table(pt[has_rare_missense == T, padjust < .5])
       # table(pt[has_rare_missense == T, padjust < 1])
       # View(pt[has_rare_missense == T])
       
       g <- ggplot(pt, aes(get(has_col), l2fc)) + geom_boxplot() + theme_bw(base_size = 14) + labs(x = has_col)
       plot(g)
       g2 <- ggplot(pt, aes(get(has_rare_col), l2fc)) + geom_boxplot() + theme_bw(base_size = 14) + labs(x = has_rare_col)
       plot(g2)
}
)




