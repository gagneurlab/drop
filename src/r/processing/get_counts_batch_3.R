suppressPackageStartupMessages(source("src/r/config.R"))
exons_op = readRDS("resources/exons_op.Rds")
bf <- substr(scan("resources/201711_nadel_rna_seq_samples.txt", what = character()), 3, 100)
samples <- vapply(strsplit(bf, "/"), "[", "", 1)

bamfiles <- BamFileList(file.path(RAWDIR, "helmholtz", bf), yieldSize=2000000)
names(bamfiles) <- samples

my_bpparam = register_bplapply_for_clustering(slurm = F, workers = 40, threads = 2, memory = 8000, jobname = "count_rna")

starttime= Sys.time()

library(GenomicAlignments)
se = summarizeOverlaps(
    exons_op,
    bamfiles,
    mode = 'IntersectionStrict',
    singleEnd = FALSE,
    ignore.strand = FALSE,    # FALSE: this batch was done strand specifically
    inter.feature = TRUE,     # TRUE, reads mapping to multiple features are dropped
    fragments = FALSE, 
    BPPARAM = my_bpparam    # SerialParam()
)
message('Processed all fibros in: ', format(Sys.time()- starttime))

saveRDS(se, "/s/project/genetic_diagnosis/processed_data/Rds/se_batch3.Rds")
message('counting done')


# Add gene info to exon count
se_ex = assay(se) %>% as.data.table()
se_ex[, exon_id := 1:.N]
exons_gene_dt <- readRDS("./resources/exons_gene_dt.Rds")
se_ex = merge(se_ex, exons_gene_dt, by = "exon_id")


# Aggregate exon data by gene and turn into a matrix
se_gene = se_ex[, lapply(.SD, sum), by = gene_name, .SDcols = colnames(assay(se))]
se_matrix = as.matrix(se_gene[, colnames(assay(se)), with = F])
row.names(se_matrix) = se_gene$gene_name
dim(se_matrix)   # 54,356 rows
se_matrix = se_matrix[rowSums2(se_matrix) > 0, ]
dim(se_matrix)   # ~ 30,000 rows

write.table(se_matrix, file.path(PROC_DATA, "rna_batch3_strand_specific.txt"), quote = F, col.names = T, row.names = T)
message('matrix done')
