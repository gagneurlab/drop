suppressPackageStartupMessages(source("src/r/config.R"))
exons_en = readRDS("resources/exons_en.Rds")

old_samples = scan("resources/kremer_bader_rna_ids.txt", what = character())

RAWDATADIR <- "/s/project/mitoMultiOmics/raw_data/helmholtz"

bamfiles = BamFileList(
    sapply(old_samples, get_helmholtz_file, 'rna','BAM.STAR'), # function to get BAM file path
    obeyQname=FALSE, 
    yieldSize= 1e6
)
names(bamfiles) = old_samples

my_bpparam = register_bplapply_for_clustering(slurm = F, workers = 40, threads = 2, memory = 8000, jobname = "count_rna")

starttime= Sys.time()

library(GenomicAlignments)
se = summarizeOverlaps(
    exons_en,
    bamfiles,
    mode = 'IntersectionStrict',
    singleEnd = FALSE,
    ignore.strand = FALSE,    # FALSE: this batch was done strand specifically
    inter.feature = TRUE,     # TRUE, reads mapping to multiple features are dropped
    fragments = FALSE, 
    BPPARAM = my_bpparam    # SerialParam()
)
message('Processed all fibros in: ', format(Sys.time()- starttime))

saveRDS(se, file.path(PROC_DATA, "se_batch0.Rds"))
message('counting done')


# Add gene info to exon count
se_ex = assay(se) %>% as.data.table()
se_ex[, exon_id := 1:.N]
exons_gene_dt <- readRDS("./resources/exons_gene_dt.Rds")
se_ex = merge(se_ex, exons_gene_dt, by = "exon_id")


# Aggregate exon data by gene and turn into a matrix
se_gene = se_ex[, lapply(.SD, sum, na.rm = T), by = gene_name, .SDcols = colnames(assay(se))]
se_matrix = as.matrix(se_gene[, colnames(assay(se)), with = F])
row.names(se_matrix) = se_gene$gene_name
dim(se_matrix)   # 54,356 rows
se_matrix = se_matrix[rowSums2(se_matrix) > 0, ]
dim(se_matrix)   # ~ 30,000 rows

write.table(se_matrix, file.path(PROC_DATA, "rna_batch0_new_annotation.txt"), quote = F, col.names = T, row.names = T)
message('matrix done')
