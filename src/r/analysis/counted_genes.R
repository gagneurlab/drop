library(magrittr)
library(data.table)
gtf_map <- readRDS("resources/GENCODEv19Mapping.RDS")
gtf_file <- readRDS("resources/gencode.v19_with_gene_name.Rds")
gtf_file <- readRDS("resources/exons_gene_dt.Rds")
dim(gtf_file)
genes_robert <- c("NDUFAF5", "NDUFAF6",
                  "NDUFA2", "NDUFA7", "NDUFA3", "NDUFA11", "NDUFA13",
                  "DNAJC30", "GCDH", "MOCS1",
                  "MRPL12", "MRPL30", "MRPL38", "MRPL45",
                  "MRPS17", "MRPS21",
                  "MSRB3", "MTG1", "RARS2", "SARS2",
                  "SDHAF2", "SLC25A10", "SLC25A26",
                  "TIMM9", "TIMM10B", "TIMM13", "TIMM23",
                  "TK2", "TOMM5", "TSFM", "ACACA", "ACAD11", "FAHD1", "GATC")

robert_ids <- gtf_file[gene_name %in% genes_robert, gene_id] %>% unique


# ods <- readRDS("/s/project/genetic_diagnosis/processed_results/ods_batches2_3_4_ss.Rds")
gene_counts <- readRDS("/s/project/genetic_diagnosis/processed_data/Rds/batches2_3_4_counts_ss.Rds")
dim(gene_counts)
gene_counts_nozero <- gene_counts[rowSums(gene_counts) > 0, ]
dim(gene_counts_nozero)

genes_robert %in% gtf_file$gene_name
genes_robert[!genes_robert %in% gtf_file$gene_name]

robert_ids %in% row.names(gene_counts)
robert_ids %in% row.names(gene_counts_nozero)

robert_ids %in% row.names(ods)

genes_robert %in% row.names(ods)
