## Parameters to change

root_dir <- '/s/project/drop-analysis/kremer_checks/run_drop_26267' # adjust
root_dir <- './'
drop_group <- 'kremer' # adjust

## OUTRIDER 
counts_out <- readRDS(file.path(root_dir, 'processed_data/aberrant_expression/v29/outrider', drop_group, 'total_counts.Rds'))
granges_out <- readRDS(file.path(root_dir, 'processed_data/preprocess/v29/count_ranges.Rds')) 
res_out <- fread(file.path(root_dir, 'processed_results/aberrant_expression/v29/outrider', drop_group, 'OUTRIDER_results.tsv'))
ods <- readRDS(file.path(root_dir, 'processed_results/aberrant_expression/v29/outrider', drop_group, 'ods.Rds'))

# FRASER 
fds_raw <- loadFraserDataSet(file = file.path(root_dir, 'processed_data/aberrant_splicing/datasets/savedObjects', paste0('raw-', drop_group), 'fds-object.RDS'))
fds <- loadFraserDataSet(file = file.path(root_dir, 'processed_results/aberrant_splicing/datasets/savedObjects', paste0(drop_group, '--v29'), 'fds-object.RDS'))
res_fra <- fread(file.path(root_dir, 'processed_results/aberrant_splicing/results/v29/fraser', drop_group, 'results.tsv'))

# MAE
mae <- fread(file.path(root_dir, 'processed_results/mae', drop_group, 'MAE_results_v29.tsv'))
qc_mat <- readRDS(file.path(root_dir, 'processed_results/mae', drop_group, 'dna_rna_qc_matrix.Rds'))



### Checks

# 1. Number of samples & genes. Use gencode29. Count unstranded.
dim(counts_out)
stopifnot(dim(counts_out) == c(length(granges_out), 119))

# 2. At least 10K expressed genes
stopifnot(nrow(ods) > 10e3)

# 3. TIMMDC1 is expression outlier in both samples & one of the pvalues
res_out[hgncSymbol == 'TIMMDC1'] # Should give 2 outliers, MUC1365 & MUC1344
stopifnot(identical(sort(res_out[hgncSymbol == 'TIMMDC1', sampleID]), c('MUC1344', 'MUC1365')))
stopifnot(assays(ods)$pValue["ENSG00000113845.9_2", "MUC1344"] == 
            res_out[geneID == 'ENSG00000113845.9_2' & sampleID == "MUC1344", pValue])

# 4. MGST1 outlier
res_out[hgncSymbol == 'MGST1'] # Should be MUC1396
stopifnot(identical(sort(res_out[hgncSymbol == 'MGST1', sampleID]), c('MUC1396')))

# 5. Check counts in one TIMMDC1 sample
stopifnot(res_out[hgncSymbol == 'TIMMDC1' & sampleID == 'MUC1344', rawcounts] == 154)

# 6. Check counts in one MGST1 sample
stopifnot(res_out[hgncSymbol == 'MGST1' & sampleID == 'MUC1396', rawcounts] == 366)


# 7. FRASER objects Dimensions
stopifnot(nrow(fds_raw) > 2e6) # Raw FRASER has > 2M rows
stopifnot(nrow(fds) > 1e5 & nrow(fds) < 5e5) # filtered FRASER has between 100-500K rows

# 8. TIMMDC1 outlier
stopifnot(identical(sort(res_fra[hgncSymbol == 'TIMMDC1', sampleID]), c('MUC1344', 'MUC1365')))

# 9. CLPP
stopifnot(identical(sort(res_fra[hgncSymbol == 'CLPP', sampleID]), c('MUC1350')))

# 10. Split counts
stopifnot(res_fra[hgncSymbol == 'CLPP' & sampleID == 'MUC1350', counts] == 40) # if this is the junction chr19 6366375 6368548  2174

# 11. MAE
stopifnot(grepl('MUC1404', mae[gene_name == 'ALDH18A1', ID]))

# 12. MAE counts
stopifnot(mae[gene_name == 'ALDH18A1' & ID == "65990--MUC1404", altCount] == 159)

# 13. DNA-RNA QC
limit <- .75
stopifnot(all(qc_mat[lower.tri(qc_mat)] < limit))
stopifnot(all(diag(qc_mat) > limit))

