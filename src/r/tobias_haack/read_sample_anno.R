source('src/r/config.R')
library(readxl)

sat <- read_xlsx("/s/project/mitoMultiOmics/raw_data/sample_info/mitOmics_samples_180606_TOH.xlsx", skip = 7) %>% as.data.table
sat <- sat[2:.N]
sat <- sat[n %in% 1:56]
setnames(sat, old = c("Con-", "No to"), new = c("Conclusion", "No_to_sequence"))
names(sat) <- gsub(" ", "_", names(sat))

sat[grepl("G$", ID_Links) | grepl("G$", Foreign_Id), ASSAY := "WGS"]   # G at the end of the ID means Genome
sat[ID_Links == "80500", ASSAY := "WES"]    # Only Exome so far
sat[grepl("R$|T$", ID_Links) | grepl("R$|T$", Foreign_Id), ASSAY := "RNASeq"]   # R/T at the end of the ID means RNA/Transcriptome
table(sat$ASSAY)
sat[, IS_RNA_SEQ_STRANDED := as.logical(IS_RNA_SEQ_STRANDED)]
sat[ASSAY == "RNASeq", table(IS_RNA_SEQ_STRANDED)]

write.table(sat, "/s/project/mitoMultiOmics/raw_data/sample_info/201812_th_sample_anno.tsv", sep = "\t", row.names = F, quote = F)
