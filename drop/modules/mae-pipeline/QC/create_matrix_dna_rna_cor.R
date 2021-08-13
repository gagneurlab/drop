#'---
#' title: Create QC matrix
#' author: vyepez
#' wb:
#'  log:
#'   - snakemake: '`sm str(tmp_dir / "MAE" / "{dataset}" / "QC_matrix.Rds")`'
#'  params:
#'    - rnaIds: '`sm lambda w: sa.getIDsByGroup(w.dataset, assay="RNA")`'
#'  input: 
#'    - mae_res: '`sm lambda w: expand(cfg.getProcessedDataDir() +
#'                "/mae/RNA_GT/{rna}.Rds", rna=sa.getIDsByGroup(w.dataset, assay="RNA"))`'
#'  output:
#'    - mat_qc: '`sm cfg.getProcessedResultsDir() + 
#'               "/mae/{dataset}/dna_rna_qc_matrix.Rds"`'
#'  threads: 20
#'  type: script
#'---

saveRDS(snakemake, snakemake@log$snakemake)

suppressPackageStartupMessages({
  library(VariantAnnotation)
  library(GenomicRanges)
  library(magrittr)
  library(BiocParallel)
  library(data.table)
})

register(MulticoreParam(snakemake@threads))
sa <- fread(snakemake@config$sampleAnnotation)

# Read the test vcf as GRanges
gr_test <- readVcf(snakemake@config$mae$qcVcf) %>% granges()
mcols(gr_test)$GT <- "0/0"

# Obtain the rna and vcf files
rna_samples <- snakemake@params$rnaIds
mae_res <- snakemake@input$mae_res

vcf_cols <- sa[RNA_ID %in% rna_samples, .(DNA_ID, DNA_VCF_FILE)] %>% unique %>% na.omit()
wes_samples <- vcf_cols$DNA_ID
vcf_files <- vcf_cols$DNA_VCF_FILE


N <- length(vcf_files)
lp <- bplapply(1:N, function(i){
  
  # Read sample vcf file
  sample <- wes_samples[i]
  param <-  ScanVcfParam(fixed=NA, info='NT', geno='GT', samples=sample, trimEmpty=TRUE) 
  vcf_sample <- readVcf(vcf_files[i], param = param, row.names = FALSE)
  # Get GRanges and add Genotype
  gr_sample <- granges(vcf_sample)
  
  if(!is.null(geno(vcf_sample)$GT)){
    gt <- geno(vcf_sample)$GT
    gt <- gsub('0|0', '0/0', gt, fixed = TRUE)
    gt <- gsub('0|1', '0/1', gt, fixed = TRUE)
    gt <- gsub('1|0', '0/1', gt, fixed = TRUE)
    gt <- gsub('1|1', '1/1', gt, fixed = TRUE)
  } else if(!is.null(info(vcf_sample)$NT)){
    gt <- info(vcf_sample)$NT
    gt <- gsub('ref', '0/0', gt)
    gt <- gsub('het', '0/1', gt)
    gt <- gsub('hom', '1/1', gt)
 }
  
  mcols(gr_sample)$GT <- gt
  
  # Find overlaps between test and sample
  gr_res <- copy(gr_test)
  seqlevelsStyle(gr_res) <- seqlevelsStyle(seqlevelsInUse(gr_sample))  # Make chr style the same
  ov <- findOverlaps(gr_res, gr_sample, type = 'equal')
  mcols(gr_res)[from(ov),]$GT <- mcols(gr_sample)[to(ov),]$GT
  
  # Find simmilarity between DNA sample and RNA sample
  x <- vapply(mae_res, function(m){
    gr_rna <- readRDS(m)
    seqlevelsStyle(gr_rna) <- seqlevelsStyle(gr_res)
    ov <- findOverlaps(gr_res, gr_rna, type = 'equal')
    gt_dna <- gr_res[from(ov)]$GT
    gt_rna <- gr_rna[to(ov)]$RNA_GT
    sum(gt_dna == gt_rna) / length(gt_dna)
  }, 1.0)
  return(x)
})

# Create a matrix
mat <- do.call(rbind, lp)
row.names(mat) <- wes_samples
colnames(mat) <- rna_samples


saveRDS(mat, snakemake@output$mat_qc)

