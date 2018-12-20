source("src/r/config.R")
source("src/r/tobias_haack/read_sample_anno.R")
library(OUTRIDER)
library(AnnotationDbi)

# Read annotation file
gencode_txdb <- loadDb("../genetic_diagnosis/resources/gencode.v19.genes.patched_contigs.Db")

# Read all batches
batches_nss <- c(
    "se_batch0.Rds",
    "se_batch1_non_strand_specific.Rds")

batches_ss <- c(
    "se_batch2_strand_specific.Rds",
    "se_batch3.Rds",
    "se_batch_arc.Rds"
)

# Strand specific samples
ss_samples <- SAMPLE_ANNOTATION[!is.na(RNA_ID) & 
                                    TISSUE == "FIBROBLAST" & 
                                    GROWTH_MEDIUM == "GLU" &
                                    is.na(TRANSDUCED_GENE) &
                                    IS_RNA_SEQ_STRANDED == T, RNA_ID]
ss_samples <- c(ss_samples, sat[IS_RNA_SEQ_STRANDED == T, ID_Links])  # Add Tobias's samples
ss_samples %>% length()

# Non strand specific samples
nss_samples <- SAMPLE_ANNOTATION[!is.na(RNA_ID) & 
                                    TISSUE == "FIBROBLAST" & 
                                    GROWTH_MEDIUM == "GLU" &
                                    is.na(TRANSDUCED_GENE) &
                                    IS_RNA_SEQ_STRANDED == F, RNA_ID]
nss_samples <- c(nss_samples, sat[IS_RNA_SEQ_STRANDED == F, ID_Links])  # Add Tobias's samples
nss_samples %>% length()


DIR_rna_samples <- "/s/project/genetic_diagnosis/processed_data/Rds/samples"
DIR_Rds <- "/s/project/genetic_diagnosis/processed_data/Rds"

rna_ss_files <- c(intersect(list.files(DIR_rna_samples, full.names = T), paste0(DIR_rna_samples, "/", ss_samples, ".Rds")),
                   paste0(file.path(DIR_Rds, batches_ss)))

rna_nss_files <- c(intersect(list.files(DIR_rna_samples, full.names = T), paste0(DIR_rna_samples, "/", nss_samples, ".Rds")),
                  paste0(file.path(DIR_Rds, batches_nss)))


gene_counts_from_se <- function(se){
    exons_gene_dt <- readRDS("../genetic_diagnosis/resources/exons_gene_dt.Rds")
    
    # Add gene level info
    exon_counts_dt <- as.data.table(se)
    exon_counts_dt[, exon_id := 1:.N]
    exon_counts_dt = merge(exon_counts_dt, exons_gene_dt, by = "exon_id")
    
    # Aggregate exon data by gene and turn into a matrix
    gene_counts_dt = exon_counts_dt[, lapply(.SD, sum), by = gene_id, .SDcols = colnames(se)]
    gene_counts = as.matrix(gene_counts_dt[, colnames(se), with = F])
    row.names(gene_counts) = gene_counts_dt[, gene_id]
    gene_counts
}


# Read RDS files, take the counts and cbind them
ss_counts <- do.call(cbind, lapply(lapply(rna_ss_files, readRDS), function(x) assay(x, "counts")))

bad_samples <- setdiff(SAMPLE_ANNOTATION[, unique(RNA_ID)], ss_samples)
ss_counts <- ss_counts[, ! colnames(ss_counts) %in% bad_samples]

ss_counts_gene <- gene_counts_from_se(ss_counts)
ss_counts_gene <- ss_counts_gene[rowSums(ss_counts_gene) > 0, ]
dim(ss_counts_gene)

saveRDS(ss_counts_gene, "/s/project/genetic_diagnosis/processed_data/Rds/batches2_3_4_th_counts_ss.Rds")
s1 <- readRDS("/s/project/genetic_diagnosis/processed_data/Rds/batches2_3_4_counts_ss.Rds")
dim(s1)

# Read RDS files, take the counts and cbind them
nss_counts <- do.call(cbind, lapply(lapply(rna_nss_files, readRDS), function(x) assay(x, "counts")))

bad_samples <- setdiff(SAMPLE_ANNOTATION[, unique(RNA_ID)], nss_samples)
nss_counts <- nss_counts[, ! colnames(nss_counts) %in% bad_samples]

nss_counts_gene <- gene_counts_from_se(nss_counts)
nss_counts_gene <- nss_counts_gene[rowSums(nss_counts_gene) > 0, ]

saveRDS(nss_counts_gene, "/s/project/genetic_diagnosis/processed_data/Rds/batches0_1_th_counts_ss.Rds")
# nss_counts_gene <- readRDS("/s/project/genetic_diagnosis/processed_data/Rds/batches0_1_th_counts_ss.Rds")

dim(nss_counts_gene)
colnames(nss_counts_gene)

# create outrider object
library(OUTRIDER)
ods <- OutriderDataSet(countData = ss_counts_gene)
ods <- OutriderDataSet(countData = nss_counts_gene)
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
g
ggsave("/s/project/genetic_diagnosis/processed_results/hist_FPKM_3batches_ss.png", g)
ggsave("/s/project/genetic_diagnosis/processed_results/hist_FPKM_2batches_nss.png", g)

ods <- filterExpression(ods, gtfFile=gencode_txdb, filter=TRUE)
dim(ods)


# Add genes metainfo
gene_name_as_row_name <- function(ods){
    genes_dt <- readRDS("./resources/gencode.v19_with_gene_name.Rds")
    rowData(ods)$geneID = row.names(ods)
    rowData(ods) = left_join(as.data.table(rowData(ods)), genes_dt[,.(gene_id, gene_name, gene_type)], by = c("geneID" = "gene_id"))
    rownames(ods) = rowData(ods)$gene_name
    ods
}

ods_ss <- gene_name_as_row_name(ods_ss)
ods_nss <- gene_name_as_row_name(ods_nss)
ods_blood <- readRDS("/s/project/genetic_diagnosis/processed_results/ods_blood.Rds")
ods_blood <- gene_name_as_row_name(ods_blood)
saveRDS(ods_blood, "/s/project/genetic_diagnosis/processed_results/ods_blood.Rds")


# run full outrider
ods <- estimateSizeFactors(ods)
# saveRDS(counts(ods, normalized = T), "/s/project/genetic_diagnosis/processed_data/Rds/batches2_3_4_counts_ss_fpkmfiltered_normsizefactors.Rds")
saveRDS(ods, "~/Downloads/ods.Rds")
ods <- readRDS("~/Downloads/ods.Rds")

pars <- c(seq(5, min(40, ncol(ods), nrow(ods)), 2), 50, 70)
ods <- findEncodingDim(ods, lnorm = T, BPPARAM = MulticoreParam(30), params = pars)
ods <- OUTRIDER(ods, BPPARAM = MulticoreParam(30))

# ods <- autoCorrect(ods, q = 60)  # Felix recommended, q = Ngenes / 4
# ods <- fit(ods)
# ods <- computePvalues(ods)
# ods <- computeZscores(ods)

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

# saveRDS(ods, "/s/project/genetic_diagnosis/processed_results/ods_4batches.Rds")
saveRDS(ods_ss, "/s/project/genetic_diagnosis/processed_results/ods_batches2_3_4_th_ss.Rds")
saveRDS(ods_nss, "/s/project/genetic_diagnosis/processed_results/ods_batches0_1_th_nss.Rds")

ods_ss <- readRDS("/s/project/genetic_diagnosis/processed_results/ods_batches2_3_4_th_ss.Rds")

res_ss <- results(ods_ss)
res_ss[, IS_RNA_SEQ_STRANDED := T]


ods_nss <- readRDS("/s/project/genetic_diagnosis/processed_results/ods_batches0_1_th_nss.Rds")
res_nss <- results(ods_nss)
res_nss[, IS_RNA_SEQ_STRANDED := F]

res <- rbind(res_ss, res_nss)

res[, LAB := "PROKISCH"]
res[sampleID %in%  sat[, ID_Links], LAB := "HAACK"]




# do some global plotting
ods <- plotCountCorHeatmap(ods, normalized=FALSE, rowCoFactor="batch")
# ods <- ods[,!colData(ods)$clusterNumber %in% c("3", "4")]

ods <- plotCountCorHeatmap(ods, normalized=TRUE, rowCoFactor="batch")

png("/s/project/genetic_diagnosis/processed_results/barplot_aberrant_genes_per_sample.png", width=900, height=700, res=120)
plotAberrantPerSample(ods)
dev.off()

# do it if you have time and a big memory 
# plotQQ(ods, global=TRUE)

