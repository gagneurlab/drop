#'---
#' title: Variants associated to OUTRIDER
#' author: vyepez
#' wb:
#'  input:
#'   - var_dt: '`sm config["PROC_RESULTS"] + "/process_vcf/unique_vars_gene/unique_variant_dt.Rds"`'
#'   - res_fib: '`sm config["PROC_RESULTS"] + "/v29_overlap/outrider/fib/OUTRIDER_results_all.Rds"`'
#'  output:
#'   - variants_outrider: '`sm config["PROC_RESULTS"] + "/process_vcf/vars_outrider/outrider_all_vars.Rds"`'
#'   - mae_outrider: '`sm config["PROC_RESULTS"] + "/processed_results/mae/mae_results_outrider.Rds"`'
#'  threads: 10
#'  type: script
#'---

#+ echo=F
saveRDS(snakemake, "tmp/outrider_vars.snakemake")
# snakemake <- readRDS("tmp/outrider_vars.snakemake")

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

#' ## Read Data
# Read OUTRIDER results and remove samples without EXOME_IDs
sa <- fread(snakemake@config$SAMPLE_ANNOTATION)
res_fib <- readRDS(snakemake@input$res_fib)
res_fib <- left_join(res_fib, sa[GROWTH_MEDIUM == 'GLU' & is.na(TRANSDUCED_GENE), .SD[1], by = 'EXOME_ID'][,.(RNA_ID, EXOME_ID, KNOWN_MUTATION)], 
                     by = c("sampleID" = "RNA_ID")) %>% as.data.table
res_fib <- res_fib[! sampleID %in% c("102630R", "102631R", "76634", "MUC1353", "MUC1356", "MUC1366", "MUC1406", "MUC1435")]
res_fib[, geneID := toupper(geneID)]
res_fib = res_fib[!is.na(EXOME_ID)]

# Read variants table
vt <- readRDS(snakemake@input$var_dt)
# vt <- readRDS("/s/project/genetic_diagnosis/processed_results/process_vcf/unique_vars_gene/unique_variant_dt.Rds")
# vt <- readRDS("/s/project/genetic_diagnosis/processed_results/process_vcf/variant_dt.Rds")

vt[, c("var_id", "rare_mito", "quality") := NULL]
exomes_common <- intersect(unique(res_fib$EXOME_ID), unique(vt$EXOME_ID))
vt = vt[EXOME_ID %in% exomes_common]
vt[, mstype := as.character(mstype)]

# Merge outrider results with variants. Remove genes with no variants.
pt <- left_join(res_fib, vt, by = c("EXOME_ID", "geneID")) %>% as.data.table
# pt <- pt[gene_type == 'protein_coding']
pt <- pt[!is.na(chr)]

#+echo=F
pt[mstype == 'ncrna-exon', mstype := 'ncrna_exon']
pt[mstype == 'frame-shift', mstype := 'frameshift']
pt[mstype == '3utr', mstype := 'utr3']
pt[mstype == '5utr', mstype := 'utr5']

#'
group_var_types <- c(splice_acceptor = "splice", splice_donor = "splice", splice = "splice",
                     frameshift = "frameshift", 
                     utr3 = "non_coding", utr5 = "non_coding", downstream = "non_coding", upstream = "non_coding", intron = "non_coding", ncrna_exon = "non_coding", mature_miRNA = "non_coding",
                     coding = "coding", del = "coding", ins = "coding", missense = "coding",
                     stop = "stop", unstop = "stop", unstart = "stop",
                     synonymous = "synonymous", stop_retain = "synonymous")

pt[, var_type := group_var_types[mstype]]
pt[, genotype := 'heterozygous']
pt[homozygous == TRUE, genotype := 'homozygous']

#' ## Asign allele frequency categories
asign_allele_freq_cat <- function(DT, breaks, column_name = 'MAX_AF', new_column_name = 'AF_CAT'){
  DT[, c(new_column_name) := NULL]   # delete the column if already present
  dt <- DT[, ..column_name]
  N <- length(breaks) - 1
  ctg <- character(N)
  for(i in 1:N){
    ctg[i] <- paste(breaks[i], breaks[i+1], sep = "-")
    dt[breaks[i] < get(column_name) & get(column_name) <= breaks[i+1], xy := ctg[i]]
  }
  dt[, xy := factor(xy, levels = ctg)]
  setnames(dt, 'xy', new_column_name)
  DT <- data.table(DT, dt[, ..new_column_name])
  return(DT)
}
pt <- asign_allele_freq_cat(pt, c(0, .001, .01, 1))

saveRDS(pt, snakemake@output$variants_outrider)



#' MAE
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

mae_res_all <- asign_allele_freq_cat(mae_res_all, c(0, .001, .01, 1))

mae_res_all[, allele_dif := 1 - 2*alt_freq]

rf <- res_fib[aberrant == T]
rf[, aux := paste(EXOME_ID, sampleID, geneID, sep = "-")]
mae_res_all[, aux := paste(sample, hgncid, sep = "-")]
mae_res_all[, aberrant := aux %in% rf$aux, by = 1:nrow(mae_res_all)]

saveRDS(mae_res_all, snakemake@output$mae_outrider)

