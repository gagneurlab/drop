# bloods <- SAMPLE_ANNOTATION[TISSUE == 'BLOOD', RNA_ID]
# bloods <- bloods[bloods != "70477"]
# bloods_files <- paste0('./', bloods, '/RNAout/paired-endout/STAR.bam')
# write.table(bloods_files, 'resources/201807_blood_rna_seq_samples.txt', sep = '\n', quote = F, row.names = F, col.names = F)


suppressPackageStartupMessages(source("src/r/config.R"))
exons_en = readRDS("resources/exons_en.Rds")  # take normal exon annotation
bf <- substr(scan("resources/201807_blood_rna_seq_samples.txt", what = character()), 3, 100)
samples <- vapply(strsplit(bf, "/"), "[", "", 1)

bamfiles <- BamFileList(file.path(RAWDIR, "helmholtz", bf), yieldSize=2000000)
names(bamfiles) <- samples

my_bpparam = register_bplapply_for_clustering(slurm = F, workers = 30, threads = 2, memory = 8000, jobname = "count_rna2")

starttime= Sys.time()

library(GenomicAlignments)
se = summarizeOverlaps(
    exons_en,
    bamfiles,
    mode = 'IntersectionStrict',
    singleEnd = FALSE,
    ignore.strand = TRUE,    # TRUE: these blood samples were done not strand specifically
    inter.feature = TRUE,     # TRUE, reads mapping to multiple features are dropped
    fragments = FALSE, 
    BPPARAM = my_bpparam    # SerialParam()
)
message('Processed all fibros in: ', format(Sys.time()- starttime))

saveRDS(se, "/s/project/genetic_diagnosis/processed_data/Rds/se_blood5.Rds")
message('counting done')


# Read all batches
# batches <- c("se_blood.Rds", "se_blood2.Rds", "se_blood3.Rds", "se_blood4.Rds", "se_blood5.Rds")

# rcounts <- lapply(batches, function(x) readRDS(file.path("/s/project/genetic_diagnosis/processed_data/Rds", x)))
# lapply(rcounts, dim)
# lapply(rcounts, colnames)

# se_blood_all = rcounts[[1]]
# for(i in 2:length(rcounts)){
#   se_blood_all = cbind(se_blood_all, rcounts[[i]])
# }
# class(se_blood_all)
# dim(se_blood_all)

se_blood <- readRDS("/s/project/genetic_diagnosis/processed_data/Rds/se_blood_all.Rds")

# Add gene info to exon count
se_ex = assay(se_blood) %>% as.data.table()
se_ex[, exon_id := 1:.N]
exons_gene_dt <- readRDS("./resources/exons_gene_dt.Rds")
se_ex = merge(se_ex, exons_gene_dt, by = "exon_id")


# Aggregate exon data by gene and turn into a matrix
se_gene = se_ex[, lapply(.SD, sum), by = gene_id, .SDcols = colnames(assay(se_blood))]
se_matrix = as.matrix(se_gene[, colnames(assay(se_blood)), with = F])
row.names(se_matrix) = se_gene$gene_id
dim(se_matrix)   # 54,356 rows
se_matrix = se_matrix[rowSums2(se_matrix) > 0, ]
dim(se_matrix)   # ~ 38,150 rows

write.table(se_matrix, file.path(PROC_DATA, "rna_blood_all.txt"), quote = F, col.names = T, row.names = T)
message('matrix done')
se_matrix <- read.delim(file.path(PROC_DATA, "rna_blood_all.txt"), sep = " ", check.names = F)


#### OUTRIDER on blood samples

# create outrider object
library(OUTRIDER)
ods_blood <- OutriderDataSet(countData = se_matrix)
colData(ods_blood)$sampleID <- colnames(ods_blood)
colData(ods_blood)$TISSUE <- "BLOOD"

# filter not expressed genes
gencode_txdb <- loadDb("../genetic_diagnosis/resources/gencode.v19.genes.patched_contigs.Db")
ods_blood <- filterExpression(ods_blood, gtfFile=gencode_txdb, filter=FALSE)
g <- plotFPKM(ods_blood) + theme_bw(base_size = 14)
ggsave("/s/project/genetic_diagnosis/processed_results/hist_blood.png", g)

ods_blood <- filterExpression(ods_blood, gtfFile=gencode_txdb, filter=TRUE)

# Add genes metainfo
rowData(ods_blood)$geneID = row.names(ods_blood)
genes_dt <- fread("/s/project/genetic_diagnosis/resource/genes_pc.txt")
rowData(ods_blood) = merge(rowData(ods_blood), genes_dt[,.(gene_id, gene_name, gene_type)], by.x = "geneID", by.y = "gene_id")
# rownames(ods_blood) = rowData(ods_blood)$gene_name

# run full outrider
ods_blood <- estimateSizeFactors(ods_blood)

ods_blood <- controlForConfounders(ods_blood, q = 40)  # Felix recommended

ods_blood <- fit(ods_blood)
ods_blood <- computePvalues(ods_blood)
ods_blood <- computeZscores(ods_blood)

ods_blood <- plotCountCorHeatmap(ods_blood, normalized=FALSE)
ods_blood <- plotCountCorHeatmap(ods_blood, normalized=TRUE)

library(biomaRt)
ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
ensgID <- rownames(ods_blood)
ensg_match <- data.table(
    ensgIDv = ensgID,
    ensgID = gsub("\\.[0-9]+", "", ensgID)
)
ensg2match <- as.data.table(getBM(c("ensembl_gene_id", "hgnc_symbol"), 
                                  filters = "ensembl_gene_id", values = ensg_match$ensgID, mart = ensembl))
ensg2match
ensg_match
res_merge <- merge(ensg2match, ensg_match, by.x="ensembl_gene_id", by.y="ensgID")
res_merge <- res_merge[order(ensembl_gene_id, hgnc_symbol)][!duplicated(ensgIDv)]
setkey(res_merge, "ensgIDv")
res_merge
toreplace <- res_merge[rownames(ods_blood)]

rownames(ods_blood) <- toreplace[,hgnc_symbol]

rownames(ods_blood)[toreplace[,hgnc_symbol == "" | is.na(hgnc_symbol)]] <- toreplace[hgnc_symbol == "" | is.na(hgnc_symbol), ensgIDv]

saveRDS(ods_blood, "/s/project/genetic_diagnosis/processed_results/ods_blood.Rds")
ods_blood <- readRDS("/s/project/genetic_diagnosis/processed_results/ods_blood.Rds")

plotAberrantPerSample(ods_blood)

res_blood <- OUTRIDER::results(ods_blood)
res_blood <- add_all_gene_info(res_blood, gene_name = "geneID")
res_blood <- res_blood[, .SD[1], by = .(geneID, sampleID)]
saveRDS(res_blood, "/s/project/genetic_diagnosis/processed_results/res_blood.Rds")
write.table(res_blood, "/s/project/genetic_diagnosis/processed_results/res_blood.txt", quote = F, col.names = T, row.names = F)
res_blood = readRDS("/s/project/genetic_diagnosis/processed_results/res_blood.Rds")

res_blood[sampleID == 'MUC1682']
plotVolcano(ods_blood, sampleID="MUC1682")

plotExpressionRank(ods_blood, geneID="PTPRC")
