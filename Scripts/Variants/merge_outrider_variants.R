#'---
#' title: Variants associated to OUTRIDER
#' author: vyepez
#' wb:
#'  input:
#'   - var_dt: '`sm config["PROC_RESULTS"] + "/process_vcf/unique_vars_gene/unique_variant_dt.Rds"`'
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

vt <- readRDS(snakemake@input$var_dt)
# vt <- readRDS("/s/project/genetic_diagnosis/processed_results/process_vcf/unique_vars_gene/unique_variant_dt.Rds")
# vt <- readRDS("/s/project/genetic_diagnosis/processed_results/process_vcf/variant_dt.Rds")

vt[, c("var_id", "rare_mito", "quality") := NULL]
exomes_common <- intersect(unique(res_fib$EXOME_ID), unique(vt$EXOME_ID))
vt = vt[EXOME_ID %in% exomes_common]
vt[, mstype := as.character(mstype)]

pt <- left_join(res_fib, vt, by = c("EXOME_ID", "geneID")) %>% as.data.table
pt <- pt[gene_type == 'protein_coding']
pt <- pt[!is.na(chr)]

pt[mstype == 'ncrna-exon', mstype := 'ncrna_exon']
pt[mstype == 'frame-shift', mstype := 'frameshift']
pt[mstype == '3utr', mstype := 'utr3']
pt[mstype == '5utr', mstype := 'utr5']

# TODO: remove mature_mRNA
group_var_types <- c(splice_acceptor = "splice", splice_donor = "splice", splice = "splice",
                     frameshift = "frameshift", 
                     utr3 = "non_coding", utr5 = "non_coding", downstream = "non_coding", upstream = "non_coding", intron = "non_coding", ncrna_exon = "non_coding", mature_miRNA = "non_coding",
                     coding = "coding", del = "coding", ins = "coding", missense = "coding",
                     stop = "stop", unstop = "stop", unstart = "stop",
                     synonymous = "synonymous", stop_retain = "synonymous")

pt[, var_type := group_var_types[mstype]]

pt[MAX_AF <= .001, AF_CAT := "<=.001"]
pt[MAX_AF > .001 & MAX_AF < .5, AF_CAT := ".001 - .5"]
pt[MAX_AF >= .5, AF_CAT := ">= .5"]
pt[, AF_CAT := factor(AF_CAT, levels = c("<=.001", ".001 - .5", ">= .5"))]

saveRDS(pt, snakemake@output$done_file)

library(ggthemes)
library(ggbeeswarm)
ggplot(pt[!is.na(AF_CAT)], aes(var_type, FC)) + geom_violin(aes(fill = AF_CAT)) + geom_hline(yintercept = 1) + 
  scale_y_log10(limits = c(.1, 10)) + scale_fill_ptol() + facet_wrap(~homozygous, nrow = 2)

ggplot(pt[!is.na(AF_CAT)], aes(var_type, FC)) + geom_boxplot(aes(fill = AF_CAT)) + geom_hline(yintercept = 1) + 
  scale_y_log10(limits = c(.1, 10)) + scale_fill_ptol() + facet_wrap(~homozygous, nrow = 2)


# outliers only
ggplot(pt[!is.na(AF_CAT) & aberrant == T], aes(var_type, FC, color = AF_CAT)) + # geom_violin(aes(fill = AF_CAT)) + 
  geom_boxplot() + 
  geom_quasirandom(dodge.width = .8) + 
  geom_hline(yintercept = 1) + 
  # scale_y_log10(limits = c(.1, 10)) +
  scale_y_log10() + 
  scale_color_ptol() + facet_wrap(~homozygous, nrow = 2) + ggtitle("Outliers only")


paper_samples_dt <- unique(pt[,.(EXOME_ID, sampleID)])
paper_samples_dt[, MAE_ID := paste(EXOME_ID, sampleID, sep = "-")]

mae_res_all <- lapply(paper_samples_dt[, MAE_ID], function(id){
  mt <- readRDS(paste0("/s/project/genetic_diagnosis/processed_results/mae/samples/", id,"_res.Rds"))
  mt[, c("noccds", "sift1", "pph1", "gnomAD_NFE_AF", "gnomAD_AFR_AF","gnomAD_EAS_AF", "gnomAD_AMR_AF", "gnomAD_ASJ_AF", "gnomAD_SAS_AF", "rsid", "pubmed") := NULL]
  mt <- mt[!is.na(hgncid)]
}) %>% rbindlist()

mae_res_all[, mstype := as.character(mstype)]

mae_res_all[mstype == 'ncrna-exon', mstype := 'ncrna_exon']
mae_res_all[mstype == 'frame-shift', mstype := 'frameshift']
mae_res_all[mstype == '3utr', mstype := 'utr3']
mae_res_all[mstype == '5utr', mstype := 'utr5']

mae_res_all[, var_type := group_var_types[mstype]]
mae_res_all[MAX_AF <= .001, AF_CAT := "<=.001"]
mae_res_all[MAX_AF > .001 & MAX_AF < .5, AF_CAT := ".001 - .5"]
mae_res_all[MAX_AF >= .5, AF_CAT := ">= .5"]
mae_res_all[, AF_CAT := factor(AF_CAT, levels = c("<=.001", ".001 - .5", ">= .5"))]

mae_res_all[, allele_dif := 1 - 2*alt_freq]

rf <- res_fib[aberrant == T]
rf[, aux := paste(EXOME_ID, sampleID, geneID, sep = "-")]
mae_res_all[, aux := paste(sample, hgncid, sep = "-")]
mae_res_all[, aberrant := aux %in% rf$aux, by = 1:nrow(mae_res_all)]

saveRDS(mae_res_all, '/s/project/genetic_diagnosis/processed_results/mae/mae_results_outrider.Rds')

ggplot(mae_res_all[!is.na(AF_CAT)], aes(var_type, allele_dif)) + geom_violin(aes(fill = AF_CAT)) + 
  scale_fill_ptol() + geom_hline(yintercept = median(mae_res_all$allele_dif))

ggplot(mae_res_all[!is.na(AF_CAT)], aes(var_type, allele_dif)) + geom_boxplot(aes(fill = AF_CAT)) + 
  scale_fill_ptol() + geom_hline(yintercept = median(mae_res_all$allele_dif), color = 'gray')


# Outliers only
ggplot(mae_res_all[!is.na(AF_CAT) & aberrant == T], aes(var_type, allele_dif, color = AF_CAT)) + 
  geom_boxplot() + 
  geom_quasirandom(dodge.width = .8) + 
  geom_hline(yintercept = median(mae_res_all$allele_dif), color = 'gray') + 
  scale_color_ptol() + ggtitle("Outliers only")

