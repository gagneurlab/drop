suppressPackageStartupMessages(source("src/r/config.R"))
exons_op = readRDS("resources/exons_op.Rds")

source("src/r/tobias_haack/read_sample_anno.R")

samples <- sat[ASSAY == "RNASeq", ID_Links]
bf <- paste0(samples, "/RNAout/paired-endout/merged.rmdup.bam")

bamfiles <- BamFileList(file.path(RAWDIR, "helmholtz", bf), yieldSize=2000000)
names(bamfiles) <- samples

my_bpparam = register_bplapply_for_clustering(slurm = F, workers = 40, threads = 2, memory = 8000, jobname = "count_rna")

# Test if the bam files exist
bplapply(1:length(bamfiles), function(i){
    file.exists(bamfiles[[i]]$path)
    
}, BPPARAM = my_bpparam
)

# Count
library(GenomicAlignments)

bplapply(1:length(bamfiles), function(i){
    se = summarizeOverlaps(
        exons_op,
        bamfiles[i],
        mode = 'IntersectionStrict',
        singleEnd = FALSE,
        ignore.strand = ! sat[ID_Links == names(bamfiles[i]), IS_RNA_SEQ_STRANDED],    # FALSE means strand specific
        inter.feature = TRUE,     # TRUE, reads mapping to multiple features are dropped
        fragments = FALSE   # SerialParam()
    )
    saveRDS(se, paste0("/s/project/genetic_diagnosis/processed_data/Rds/samples/", names(bamfiles[i]), ".Rds"))
    },
    BPPARAM = my_bpparam 
)
# message('Processed all fibros in: ', format(Sys.time()- starttime))
message('counting done')


# Do it using Daniel's method
genes_gr <- readRDS('/s/genomes/human/hg19/ucsc/ucsc_translated_genes_gr.RDS')
bplapply(1:length(bamfiles), function(i){
    se = summarizeOverlaps(
        genes_gr,
        bamfiles[i],
        mode = 'IntersectionStrict',
        singleEnd = FALSE,
        ignore.strand = T,    # count them non strand specifically to merge with previous counts
        inter.feature = TRUE,     # TRUE, reads mapping to multiple features are dropped
        fragments = FALSE   # SerialParam()
    )
    saveRDS(se, paste0("/s/project/genetic_diagnosis/processed_data/Rds/samples/", names(bamfiles[i]), "nss.Rds"))
    },
    BPPARAM = my_bpparam 
)
message('counting done')


library(OUTRIDER)

paper_counts <- fread('/s/project/mitoMultiOmics/paper_nature_genetics/paper_supplement_data/ext_suppl_rna_gene_counts_raw.tsv')
paper_matrix <- as.matrix(paper_counts[,-1])
dim(paper_matrix)  # 27682 genes
row.names(paper_matrix) <- paper_counts$V1
rm(paper_counts)

# Get Tobias samples
th_counts <- do.call(cbind, lapply(file.path('/s/project/genetic_diagnosis/processed_data/Rds/samples/', paste0(samples, "nss.Rds")), readRDS))
th_counts <- assay(th_counts)

identical(row.names(th_counts), row.names(paper_matrix))

paper_th_counts <- cbind(paper_matrix, th_counts)
dim(paper_th_counts)

# ods <- OutriderDataSet(countData = paper_th_counts)
# ods <- estimateSizeFactors(ods)
# ods <- filterExpression(ods, '/s/genomes/human/hg19/ucsc/ucsc.translated.gtf')
# ods <- OUTRIDER(ods, BPPARAM = MulticoreParam(30))

ods <- readRDS("/s/project/genetic_diagnosis/processed_results/ods_kremer_th.Rds")
res <- OUTRIDER::results(ods)
res[, LAB := 'PROKISCH']
res[sampleID %in% samples, LAB := 'HAACK']
res[LAB == 'HAACK']
res[, geneID := toupper(geneID)]

source("Scripts/_functions/gene_annotation/add_gene_info_cols.R")
res <- add_all_gene_info(res, gene_name_col = "geneID", hans = T, omim = T, mitocarta = T, dis_genes = F, rcc = F)
res[LAB == 'HAACK'] %>% View
write.table(res[LAB == 'HAACK'], "/s/project/genetic_diagnosis/processed_results/res_th.tsv", row.names = F, quote = F, sep = '\t')

