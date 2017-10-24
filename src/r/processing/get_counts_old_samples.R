ssfdsold_counts <- read.delim2("/data/ouga/home/ag_gagneur/yepez/workspace/scared/Data/sup_data1_raw_rna_gene_counts.txt")
colnames(old_counts) = gsub("X", "", colnames(old_counts))

old95 = old_counts[ apply(old_counts, 1, quantile, .95) > 10, ] 
dim(old95)  # 12606 x 119

old_samples = colnames(old_counts)

# Check if they're all called 'merged.bam'
list.files(file.path(RAWDIR, "helmholtz", unmerged_samples, "/RNAout/paired-endout/"))

merged_samples = old_samples[sapply(old_samples, function(o) "merged.bam" %in% list.files(file.path(RAWDIR, "helmholtz", o, "/RNAout/paired-endout/"))
)]
unmerged_samples = setdiff(old_samples, merged_samples)

# SAMPLE_ANNOTATION[RNA_ID %in% old_samples & FIBROBLAST_ID == "NHDF"]

bamfiles <- BamFileList(c(file.path(RAWDIR, "helmholtz", merged_samples, "/RNAout/paired-endout/merged.bam"),
                          file.path(RAWDIR, "helmholtz", unmerged_samples, "/RNAout/paired-endout/STAR.bam")), yieldSize=2000000)
names(bamfiles) <- old_samples


starttime= Sys.time()
# !!! THIS CALL IS SENSITIVE TO VERSION OF bplapply/biocparallel !!!

se_old = summarizeOverlaps(
    exons_en, 
    # genes_en,
    bamfiles, 
    mode='IntersectionStrict',
    singleEnd=FALSE,
    ignore.strand= T, # TRUE,
    inter.feature=TRUE, 	# TRUE, reads mapping to multiple features are dropped
    fragments=FALSE,
    # BPPARAM=NULL) # ,  		# should singletons, reads with unmapped pairs and 
    #   other fragments be included in counting
    BPPARAM = my_bpparam    # SerialParam()
)
message('Processed all fibros in: ', format(Sys.time()- starttime))


saveRDS(se_old, "/data/ouga/home/ag_gagneur/yepez/Desktop/se_old_gene.Rds")
### Old samples
se_old = readRDS("/data/ouga/home/ag_gagneur/yepez/Desktop/se_old_gene.Rds")
head(assay(se_old))

# Add gene info to exon count
se_exo = assay(se_old) %>% as.data.table()
se_exo[, exon_id := 1:.N]
se_exo = merge(se_exo, exons_gene_dt, by = "exon_id")

old_samples = colnames(assay(se_old))

# Aggregate exon data by gene and turn into a matrix
se_gene_o = se_exo[, lapply(.SD, sum), by = gene_name, .SDcols = old_samples]
se_matrix_o = as.matrix(se_gene_o[, old_samples, with = F])
row.names(se_matrix_o) = se_gene_o$gene_name
dim(se_matrix_o)   # 54356 rows
se_matrix_o = se_matrix_o[rowSums2(se_matrix_o) > 0, ]
dim(se_matrix_o)   # 35922 rows


write.table(se_matrix_o, file.path(PROC_DATA, "rna_batch0_new_annotation.txt"), quote = F, col.names = T, row.names = T)


se_gene_o = readRDS("/data/ouga/home/ag_gagneur/yepez/Desktop/se_old.Rds")
se_geo = assay(se_gene_o) %>% as.data.table()

