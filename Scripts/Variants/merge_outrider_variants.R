#'---
#' title: Variants associated to OUTRIDER
#' author: vyepez
#' wb:
#'  input:
#'   - var_dt: '`sm config["PROC_RESULTS"] + "/process_vcf/variant_dt.Rds"`'
#'   - res_fib: '`sm config["PROC_RESULTS"] + "/v29_overlap/outrider/fib/OUTRIDER_results_all.Rds"`'
#'  output:
#'   - done_file: '`sm config["PROC_RESULTS"] + "/process_vcf/vars_outrider/outrider_all_vars.Rds"`'
#'  threads: 10
#'  type: script
#'---

#+ echo=F
saveRDS(snakemake, "tmp/outrider_vars.snakemake")
# snakemake <- readRDS("tmp/outrider_vars.snakemake")

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(BiocParallel)
})


# Read results and remove samples without exomes
sa <- fread(snakemake@config$SAMPLE_ANNOTATION)
res_fib <- readRDS(snakemake@input$res_fib)
res_fib <- left_join(res_fib, sa[GROWTH_MEDIUM == 'GLU' & is.na(TRANSDUCED_GENE), .SD[1], by = 'EXOME_ID'][,.(RNA_ID, EXOME_ID, KNOWN_MUTATION)], 
                     by = c("sampleID" = "RNA_ID")) %>% as.data.table
res_fib[, geneID := toupper(geneID)]
res_fib = res_fib[!is.na(EXOME_ID)]

vt <- readRDS("/s/project/genetic_diagnosis/processed_results/process_vcf/unique_vars_gene/unique_variant_dt.Rds")
vt[, c("var_id", "rare_mito", "quality") := NULL]
vt[, mstype := as.character(mstype)]
exomes_common <- intersect(unique(res_fib$EXOME_ID), unique(vt$EXOME_ID))
vt = vt[EXOME_ID %in% exomes_common]

pt <- left_join(res_fib[EXOME_ID %in% exomes_common], vt, by = c("EXOME_ID", "geneID")) %>% as.data.table
pt <- pt[gene_type == 'protein_coding']

register(MulticoreParam(snakemake@threads))

var_types <- c("missense", "synonymous", "splice", "unstop", "frame-shift", "unstart", "stop", "stop_retain")

# Merge variants table with all outrider results
bplapply(var_types,
  function(var_type){
       mt <- readRDS(paste0("/s/project/genetic_diagnosis/processed_results/process_vcf/unique_vars_gene/", var_type, "_unique_variant_dt.Rds"))
       exomes_common <- intersect(res_fib$EXOME_ID, mt$EXOME_ID)
       length(exomes_common)
       # Join only samples that are both on OUTRIDER and with variants
       pt <- left_join(res_fib[EXOME_ID %in% exomes_common], mt, by = c("EXOME_ID", "geneID")) %>% as.data.table
       
       has_col <- paste0("has_", var_type)
       has_rare_col <- paste0("has_rare_", var_type)
       set(pt, j = has_col, value = !is.na(pt$pos))
       set(pt, j = has_rare_col, value = !is.na(pt$pos) & pt$MAX_AF < 0.001)
       saveRDS(pt, paste0("/s/project/genetic_diagnosis/processed_results/process_vcf/vars_outrider/", var_type, ".Rds"))
    })


pt <- readRDS("/s/project/genetic_diagnosis/processed_results/process_vcf/vars_outrider/missense.Rds")
PT <- pt[, .(geneID,sampleID,pValue,padjust,zScore,l2fc,rawcounts,normcounts,meanCorrected,theta,aberrant,AberrantBySample,AberrantByGene,
             padj_rank,FC,gene_type,EXOME_ID,KNOWN_MUTATION,has_missense,has_rare_missense)]
for(var_type in var_types[var_types!='missense']){
  pt <- readRDS(paste0("/s/project/genetic_diagnosis/processed_results/process_vcf/vars_outrider/", var_type, ".Rds"))
  has_col <- paste0("has_", var_type)
  has_rare_col <- paste0("has_rare_", var_type)
  PT <- left_join(PT, pt[, c("geneID", "sampleID", has_col, has_rare_col), with = F], by = c("geneID", "sampleID")) %>% as.data.table
}

saveRDS(PT, snakemake@output$done_file)
