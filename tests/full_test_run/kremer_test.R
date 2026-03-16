#!/usr/env Rscript
#
# Rscript to test that we can recall all the known outliers and events
# in the Kremer et. al. dataset after a successful DROP run.
#

## Parameters to change
suppressPackageStartupMessages({
  library(data.table)
  library(OUTRIDER)
  library(FRASER)
  library(yaml)
})

# default groups and folders
root_dir <- './'
outrider_drop_group1 <- 'group1'
outrider_drop_group2 <- 'group2'
fraser_drop_group <- 'fraser'
mae_drop_group <- 'mae'

## Sample Annotation
sample_anno <- fread(read_yaml(file.path(root_dir, "config.yaml"))$sampleAnnotation)

## OUTRIDER 
counts_out <- readRDS(file.path(root_dir, 'processed_data/aberrant_expression/v29/outrider', outrider_drop_group2, 'total_counts.Rds'))
granges_out <- readRDS(file.path(root_dir, 'processed_data/aberrant_expression/v29/count_ranges.Rds'))
res_out <- fread(file.path(root_dir, 'processed_results/aberrant_expression/v29/outrider', outrider_drop_group2, 'OUTRIDER_results.tsv'))
ods <- readRDS(file.path(root_dir, 'processed_results/aberrant_expression/v29/outrider', outrider_drop_group2, 'ods.Rds'))
ods1 <- readRDS(file.path(root_dir, 'processed_results/aberrant_expression/v29/outrider', outrider_drop_group1, 'ods.Rds'))


# FRASER 
fds_raw <- loadFraserDataSet(file = file.path(root_dir, 'processed_data/aberrant_splicing/datasets/savedObjects', paste0('raw-', fraser_drop_group), 'fds-object.RDS'))
fds <- loadFraserDataSet(file = file.path(root_dir, 'processed_results/aberrant_splicing/datasets/savedObjects', paste0(fraser_drop_group, '--v29'), 'fds-object.RDS'))
res_fra <- fread(file.path(root_dir, 'processed_results/aberrant_splicing/results/v29/fraser', fraser_drop_group, 'results.tsv'))


# MAE
mae <- fread(file.path(root_dir, 'processed_results/mae', mae_drop_group, 'MAE_results_v29.tsv'))
qc_mat <- readRDS(file.path(root_dir, 'processed_results/mae', mae_drop_group, 'dna_rna_qc_matrix.Rds'))


### Checks
##
## OUTRIDER
##

# 1. Number of samples & genes. Use gencode29. Count unstranded.
dim(counts_out)
stopifnot(dim(counts_out) == c(length(granges_out), sum(sample_anno[, grepl(outrider_drop_group2, DROP_GROUP)])))
stopifnot(ncol(ods1) == sum(sample_anno[, grepl(outrider_drop_group1, DROP_GROUP)]))

# 2. At least 10K expressed genes
stopifnot(nrow(ods) > 10e3 && nrow(ods) < 20e3)

# 3. TIMMDC1 is expression outlier in both samples & one of the pvalues
res_out[hgncSymbol == 'TIMMDC1'] # Should give 2 outliers, MUC1365 & MUC1344
stopifnot(identical(sort(res_out[hgncSymbol == 'TIMMDC1', sampleID]), c('MUC1344', 'MUC1365')))
stopifnot(all.equal(assays(ods)$pValue["ENSG00000113845.9_2", "MUC1344"],
            res_out[geneID == 'ENSG00000113845.9_2' & sampleID == "MUC1344", pValue]))

# 4. MGST1 outlier
res_out[hgncSymbol == 'MGST1'] # Should be MUC1396
stopifnot(identical(sort(res_out[hgncSymbol == 'MGST1', sampleID]), c('MUC1396')))

# 5. Check counts in one TIMMDC1 sample
stopifnot(res_out[hgncSymbol == 'TIMMDC1' & sampleID == 'MUC1344', rawcounts] == 154)

# 6. Check counts in one MGST1 sample
stopifnot(res_out[hgncSymbol == 'MGST1' & sampleID == 'MUC1396', rawcounts] == 366)

# 6a. check result dimensions
stopifnot(all(sapply(c("MUC1344", "MUC1396", "MUC1365"), function(x) nrow(res_out[sampleID == x])) < c(10, 15, 15)))


### Checks
##
## FRASER
##

# 7. FRASER objects Dimensions
stopifnot(nrow(fds_raw) > 2e6) # Raw FRASER has > 2M rows
stopifnot(nrow(fds) > 85000 & nrow(fds) < 2e5) # filtered FRASER has between 80-200K rows for 80 samples
stopifnot(ncol(fds) == sum(sample_anno[, grepl(fraser_drop_group, DROP_GROUP)]))

# 7a. check result dimensions per sample
stopifnot(all(sapply(c("MUC1344", "MUC1350", "MUC1365"), function(x) nrow(res_fra[sampleID == x])) < c(15, 10, 25)))

# 8. TIMMDC1 outlier
stopifnot(identical(sort(res_fra[hgncSymbol == 'TIMMDC1', sampleID]), c('MUC1344', 'MUC1365')))

# 9. CLPP
stopifnot(identical(sort(res_fra[hgncSymbol == 'CLPP', sampleID]), c('MUC1350')))

# 10. Split counts
stopifnot(res_fra[hgncSymbol == 'CLPP' & sampleID == 'MUC1350', counts] == 36) # junction chr19 6366375 6368548  2174
stopifnot(res_fra[hgncSymbol == 'TIMMDC1' & sampleID == 'MUC1344', counts] == 37) # junction chr3 119234787 119236051  1265


### Checks
##
## MAE
##

# 11. MAE
stopifnot(nrow(mae[gene_name == 'ALDH18A1' & ID == '65990--MUC1404']) == 1)

# 12. MAE counts
stopifnot(mae[gene_name == 'ALDH18A1' & ID == "65990--MUC1404", altCount] == 159)

# 13. DNA-RNA QC report any mismatch with the annotation
limit <- .75
qc_matches <- which(qc_mat > limit, arr.ind=TRUE)
qc_match_anno <- sapply(seq_len(nrow(qc_matches)), function(x){
    rna_id = colnames(qc_mat)[qc_matches[x,'col']]
    dna_id = rownames(qc_mat)[qc_matches[x,'row']]
    ans = nrow(sample_anno[RNA_ID == rna_id & DNA_ID == dna_id]) == 1
    if(! ans){
        message("RNA_ID does not match DNA_ID:\t", rna_id, "\t", dna_id)
    }
    ans
})
stopifnot(qc_match_anno)


