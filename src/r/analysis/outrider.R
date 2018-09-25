source("src/r/config.R")
library(OUTRIDER)
library(AnnotationDbi)
library(cowplot)

# Read all batches
batches <- c(
    "se_batch0.Rds",
    "se_batch1_non_strand_specific.Rds",
    "se_batch2_strand_specific.Rds",
    "se_batch3.Rds",
    "se_batch_arc.Rds"
)


rcounts <- lapply(batches, function(x) readRDS(file.path("/s/project/genetic_diagnosis/processed_data/Rds", x)))
lapply(rcounts, dim)
lapply(rcounts, colnames)

allrcounts <- do.call(cbind, lapply(rcounts, function(x) assay(x, "counts")))
dim(allrcounts)

gencode_txdb <- loadDb("../genetic_diagnosis/resources/gencode.v19.genes.patched_contigs.Db")
exons_gene_dt <- readRDS("../genetic_diagnosis/resources/exons_gene_dt.Rds")

# Change from exon level to gene level
allrcounts_dt <- as.data.table(allrcounts)
allrcounts_dt[, exon_id := 1:.N]
allrcounts_dt = merge(allrcounts_dt, exons_gene_dt, by = "exon_id")

# Aggregate exon data by gene and turn into a matrix
allrcounts_gene_dt = allrcounts_dt[, lapply(.SD, sum), by = gene_id, .SDcols = colnames(allrcounts)]
allrcounts_gene = as.matrix(allrcounts_gene_dt[, colnames(allrcounts), with = F])
row.names(allrcounts_gene) = allrcounts_gene_dt[, gene_id]

saveRDS(allrcounts_gene, "/s/project/genetic_diagnosis/processed_data/Rds/all4batches_counts.Rds")
allrcounts_gene <- readRDS("/s/project/genetic_diagnosis/processed_data/Rds/all4batches_counts.Rds")
dim(allrcounts_gene)

# Remove other tissues and galactose
valid_samples <- SAMPLE_ANNOTATION[!is.na(RNA_ID) & TISSUE == "FIBROBLAST" & 
                                       GROWTH_MEDIUM == "GLU", RNA_ID]
cols <- intersect(colnames(allrcounts_gene), valid_samples)
allrcounts_gene <- allrcounts_gene[, cols]
dim(allrcounts_gene)

# create outrider object
library(OUTRIDER)
ods <- OutriderDataSet(countData = allrcounts_gene)
colData(ods)$sampleID <- colnames(ods)

# Add batches for heatmap
colData(ods)$batch <- as.character(NA)
for(i in seq_along(batches)){
    idxSampleBatch <- colnames(ods) %in% colnames(rcounts[[i]])
    colData(ods)$batch[idxSampleBatch] <- batches[i]
}

# filter not expressed genes
ods <- filterExpression(ods, gtfFile=gencode_txdb, filter=FALSE)
g <- plotFPKM(ods) + theme_bw(base_size = 14)
ggsave("/s/project/genetic_diagnosis/processed_results/hist_FPKM_4batches.png", g)

ods <- filterExpression(ods, gtfFile=gencode_txdb, filter=TRUE)

# Add genes metainfo
rowData(ods)$geneID = row.names(ods)
genes_dt <- fread("/s/project/genetic_diagnosis/resource/genes_pc.txt")
rowData(ods) = merge(rowData(ods), genes_dt[,.(gene_id, gene_name, gene_type)], by.x = "geneID", by.y = "gene_id")
rownames(ods) = rowData(ods)$gene_name

# run full outrider
ods <- estimateSizeFactors(ods)
ods <- autoCorrect(ods, q = 60)  # Felix recommended, q = Ngenes / 4

ods <- fit(ods)
ods <- computePvalues(ods)
ods <- computeZscores(ods)

library(biomaRt)
ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
ensgID <- rownames(ods)
ensg_match <- data.table(
    ensgIDv = ensgID,
    ensgID = gsub("\\.[0-9]+", "", ensgID)
)
ensg2match <- as.data.table(getBM(c("ensembl_gene_id", "hgnc_symbol"), filters = "ensembl_gene_id", values = ensg_match$ensgID, mart = ensembl))
ensg2match
ensg_match
res_merge <- merge(ensg2match, ensg_match, by.x="ensembl_gene_id", by.y="ensgID")
res_merge <- res_merge[order(ensembl_gene_id, hgnc_symbol)][!duplicated(ensgIDv)]
setkey(res_merge, "ensgIDv")
res_merge
toreplace <- res_merge[rownames(ods)]

rownames(ods) <- toreplace[,hgnc_symbol]

rownames(ods)[toreplace[,hgnc_symbol == "" | is.na(hgnc_symbol)]] <- toreplace[hgnc_symbol == "" | is.na(hgnc_symbol), ensgIDv]

saveRDS(ods, "/s/project/genetic_diagnosis/processed_results/ods_4batches.Rds")


# do some global plotting
ods <- plotCountCorHeatmap(ods, normalized=FALSE, rowCoFactor="batch")
# ods <- ods[,!colData(ods)$clusterNumber %in% c("3", "4")]

ods <- plotCountCorHeatmap(ods, normalized=TRUE, rowCoFactor="batch")

png("/s/project/genetic_diagnosis/processed_results/barplot_aberrant_genes_per_sample.png", width=900, height=700, res=120)
plotAberrantPerSample(ods)
dev.off()

# do it if you have time and a big memory 
# plotQQ(ods, global=TRUE)

