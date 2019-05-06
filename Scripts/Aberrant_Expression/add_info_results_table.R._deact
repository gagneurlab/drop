# Add extra columns to the OUTRIDER results table
source("src/r/config.R")
source("src/r/functions/gene_annotation/add_gene_info_cols.R")
library(OUTRIDER)

### TODO:
# Check id mismatches

# Read the OUTRIDER dataset object
ods <- readRDS("/s/project/genetic_diagnosis/processed_results/ods_4batches.Rds")
ods <- readRDS("/s/project/genetic_diagnosis/processed_results/ods_batches2_3_4_ss.Rds")

# Gene information table
genes_dt <- readRDS("resources/gencode.v19_with_gene_name.Rds")

# Subset sample annotation
sa <- SAMPLE_ANNOTATION[LAB == "PROKISCH" & !is.na(RNA_ID)]


# Get the results table
res <- OUTRIDER::results(ods, all = F)

res <- add_all_gene_info(res, gene_name_col = "geneID", hans = T, omim = T, mitocarta = T, dis_genes = F, rcc = F)
# res <- add_all_gene_info(res, gene_name_col = "gene_name", hans = T, omim = T, mitocarta = T, dis_genes = F, rcc = F)


# Add sample annotation
res <- left_join(res, sa[, .(RNA_ID, FIBROBLAST_ID, EXOME_ID, PEDIGREE, KNOWN_MUTATION,
                             CANDIDATE_GENE, BIOCHEMICAL_DEFECT, CLINICAL_SYMPTOMS)],
                 by = c("sampleID" = "RNA_ID")) %>% as.data.table


setorder(res, AberrantBySample)
dim(res)

# saveRDS(res, "/s/project/genetic_diagnosis/processed_results/res_4batches.Rds")
saveRDS(res, "/s/project/genetic_diagnosis/processed_results/res_batches2_3_4_ss.Rds")
saveRDS(res, "/s/project/genetic_diagnosis/processed_results/res_all_batches_th.Rds")
write.table(res, "/s/project/genetic_diagnosis/processed_results/res_all_batches_th.tsv", sep = "\t", row.names = F, quote = F)
res <- readRDS("/s/project/genetic_diagnosis/processed_results/res_4batches.Rds")
