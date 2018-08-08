# Add extra columns to the OUTRIDER results table
source("src/r/config.R")
source("src/r/functions/gene_annotation/add_gene_info_cols.R")
library(OUTRIDER)

### TODO:
# Check id mismatches


# Read the OUTRIDER dataset object
ods <- readRDS("/s/project/genetic_diagnosis/processed_results/ods_4batches.Rds")

# Gene information table
genes_dt <- fread("/s/project/genetic_diagnosis/resource/genes_pc.txt")

# Subset sample annotation
sa <- SAMPLE_ANNOTATION[LAB == "PROKISCH" & !is.na(RNA_ID)]


# Get the results table
res <- OUTRIDER::results(ods, all = F)
res[, c("mu", "disp") := NULL]

res <- add_all_gene_info(res, gene_name_col = "geneID")


# Add sample annotation
res <- left_join(res, sa[, .(RNA_ID, FIBROBLAST_ID, EXOME_ID, PEDIGREE, KNOWN_MUTATION,
                             CANDIDATE_GENE, BIOCHEMICAL_DEFECT, CLINICAL_SYMPTOMS)],
                 by = c("sampleID" = "RNA_ID")) %>% as.data.table


setorder(res, AberrantBySample)
dim(res)

saveRDS(res, "/s/project/genetic_diagnosis/processed_results/res_4batches.Rds")
