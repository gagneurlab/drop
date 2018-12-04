library(data.table)
library(magrittr)
sa = fread("../sample_annotation/Data/sample_annotation.tsv")[LAB == "PROKISCH"]
sa_protein = sa[! is.na(PROTEOME_ID)]
sa_sub = sa_protein[! grep("-T|_2w|yes", PROTEOME_ID)]

annot_proteomes <- sa_sub$PROTEOME_ID %>% unique

load("/s/project/mitoMultiOmics/processed_expression/proteome_normalize_fibro_bundle.Rdata")
dim(prot_100min_fibro_lfq)   # 29 x 5431
prot_100 <- colnames(prot_100min_fibro_lfq)
setdiff(prot_100, annot_proteomes)
sa[EXOME_ID %in% setdiff(prot_100, annot_proteomes)]

dim(prot_60min_fibro_lfq)    # 15 x 3517
length(proteome_limma_res_list)
prot_60 <- colnames(prot_60min_fibro_lfq)
setdiff(prot_60, annot_proteomes)  # NHDF only
sa[EXOME_ID %in% setdiff(prot_100, annot_proteomes)]

intersect(prot_100, prot_60)   # No repeated samples

setdiff(union(prot_100, prot_60), names(proteome_limma_res_list))  # NHDF not present in results



library(readxl)
prot_robert_dt <- read_excel("/s/project/mitoMultiOmics/raw_data/proteome/20180313_robert_proteomics/normalization/TMT1_6_row_col_norm.xlsx") %>% as.data.table()
prot_robert <- as.matrix(prot_robert_dt[, 2:ncol(prot_robert_dt)])
row.names(prot_robert) <- prot_robert_dt$ProteinID
prot_robert[1:5, 1:5]
prot_new <- colnames(prot_robert)
prot_new = prot_new[- grep("\\.", prot_new)]

setdiff(prot_new, annot_proteomes)

