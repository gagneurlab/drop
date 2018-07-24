# Add extra columns to the OUTRIDER results table
source("src/r/config.R")
library(OUTRIDER)

# Read the OUTRIDER dataset object
ods <- readRDS("/s/project/genetic_diagnosis/processed_results/ods_4batches.Rds")

# Gene information table
genes_dt <- fread("/s/project/genetic_diagnosis/resource/genes_pc.txt")

# Subset sample annotation
sa <- SAMPLE_ANNOTATION[LAB == "PROKISCH" & !is.na(RNA_ID)]

# Read OMIM table
omim_dt <- readRDS("../mitomultiomics/resource/omim_dt.Rds")
omim_dt <- omim_dt[SYMBOL != "", .(SYMBOL, GMIM, PINH, PMIM)]

# Get the results table
res <- results(ods, all = F)
res[, c("mu", "disp") := NULL]

# Add gene id and gene type (protein coding, antisense, etc)
res = merge(res, genes_dt[,.(gene_id, gene_name, gene_type)], by.x = "geneID", by.y = "gene_name")

# Add sample annotation
res <- left_join(res, sa[, .(RNA_ID, FIBROBLAST_ID, EXOME_ID, PEDIGREE, KNOWN_MUTATION,
                             CANDIDATE_GENE, BIOCHEMICAL_DEFECT, CLINICAL_SYMPTOMS)],
                 by = c("sampleID" = "RNA_ID")) %>% as.data.table

res <- left_join(res, omim_dt, 
                  by = c("geneID" = "SYMBOL")) %>% as.data.table()

setorder(res, AberrantBySample)
dim(res)
