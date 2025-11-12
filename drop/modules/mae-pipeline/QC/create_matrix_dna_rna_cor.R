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
#'  threads: 30
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
sa <- fread(snakemake@config$sampleAnnotation, 
            colClasses = c(RNA_ID = 'character', DNA_ID = 'character'))

# Read the test vcf as GRanges
gr_test <- readVcf(snakemake@config$mae$qcVcf) %>% granges()
mcols(gr_test)$GT <- "0/0"

# Obtain the rna and vcf files
rna_samples <- as.character(snakemake@params$rnaIds)
mae_res <- snakemake@input$mae_res

rows_in_group <- sapply(strsplit(sa$DROP_GROUP, ',|, '), function(d) snakemake@wildcards$dataset %in% d)
vcf_cols <- sa[rows_in_group, .(DNA_ID, DNA_VCF_FILE)] %>% unique
dna_samples <- vcf_cols$DNA_ID
vcf_files <- vcf_cols$DNA_VCF_FILE

# Read all RNA genotypes into a list
rna_gt_list <- bplapply(mae_res, readRDS)

# For every DNA, find the overlapping variants with each RNA
N <- length(vcf_files)

lp <- bplapply(1:N, function(i){
  
  gr_res <- copy(gr_test)
  
  # Read sample vcf file
  
  dna_id <- dna_samples[i] %>% as.character()
  
  param <-  ScanVcfParam(fixed=NA, info='NT', geno='GT', samples=dna_id, trimEmpty=TRUE) 
  vcf_dna <- readVcf(vcf_files[i], param = param, row.names = FALSE)
  # Get GRanges and add Genotype
  gr_dna <- granges(vcf_dna)
  
  if(!is.null(geno(vcf_dna)$GT)){
    gt <- geno(vcf_dna)$GT
    
    gt <- gsub('|', '/', gt, fixed = TRUE)
    gt <- gsub('1/0', '0/1', gt, fixed = TRUE)
  } else if(!is.null(info(vcf_dna)$NT)){
    gt <- info(vcf_dna)$NT
    gt <- gsub('ref', '0/0', gt)
    gt <- gsub('het', '0/1', gt)
    gt <- gsub('hom', '1/1', gt)
  }
  
  mcols(gr_dna)$GT <- gt
  
  # Find overlaps between test and sample
  gr_res <- copy(gr_test)
  seqlevelsStyle(gr_res) <- seqlevelsStyle(seqlevelsInUse(gr_dna))[1]
  ov <- findOverlaps(gr_res, gr_dna, type = 'equal')
  mcols(gr_res)$GT[from(ov)] <- mcols(gr_dna)$GT[to(ov)]
  
  # Find similarity between DNA sample and RNA sample
  sims <- sapply(rna_gt_list, function(gr_rna){
    gr_res <- gr_res[mcols(gr_res)$GT != "0/0"]
    seqlevelsStyle(gr_rna) <- seqlevelsStyle(gr_res)[1]
    ov <- findOverlaps(gr_res, gr_rna, type = 'equal')
    gt_dna <- gr_res[from(ov)]$GT
    gt_rna <- gr_rna[to(ov)]$RNA_GT
    x <- c(length(gt_rna), sum(gt_dna == gt_rna, na.rm = T))
    # x <- sum(gt_dna == gt_rna, na.rm = T)
    return(x)
    # mean(gt_dna == gt_rna, na.rm = T)
  })
  return(sims)
})

# Create a matrix
mat <- do.call(rbind, lp)
row.names(mat) <- dna_samples
colnames(mat) <- rna_samples

# Matrix 1: lengths
m1 <- do.call(rbind, lapply(lp, function(x) x[1, , drop = FALSE]))

# Matrix 2: DNA-RNA matches
m2 <- do.call(rbind, lapply(lp, function(x) x[2, , drop = FALSE]))


clean_matrices <- function(m){
  row.names(m) <- dna_samples
  colnames(m) <- rna_samples
  # Sort as in the original sa and take care of repeated DNAs/RNAs in the group
  m <- m[unique(sa[rows_in_group, DNA_ID]), 
        unique(sa[rows_in_group, RNA_ID]),
        drop=FALSE]
  
  return(m)
}

m1 <- clean_matrices(m1)
m2 <- clean_matrices(m2)

# Matrix 3: DNA-RNA matches prop.
mat <- m2/m1


saveRDS(mat, snakemake@output$mat_qc)
saveRDS(m1, gsub('matrix.Rds', 'matrix_length.Rds', snakemake@output$mat_qc))
saveRDS(m2, gsub('matrix.Rds', 'matrix_sum.Rds', snakemake@output$mat_qc))
