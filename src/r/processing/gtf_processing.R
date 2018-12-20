# Author: mumichae, vyepez

suppressPackageStartupMessages(source("src/r/config.R"))

## Annotations. Already created, just load it.
# v19 downloaded from https://storage.googleapis.com/gtex_analysis_v7/reference/gencode.v19.genes.v7.patched_contigs.gtf
# gtex_txdb = makeTxDbFromGFF('./resources/gencode.v19.genes.patched_contigs.gtf.gz', format='gtf')
# saveDb(gtex_txdb, "./resources/gencode.v19.genes.patched_contigs.Db")
# rename chromosomes for v19
gtex_txdb <- loadDb("./resources/gencode.v19.genes.patched_contigs.Db")
seqlevels(gtex_txdb) <- paste0("chr", seqlevels(gtex_txdb))
seqlevels(gtex_txdb) <- gsub("MT", "M", seqlevels(gtex_txdb))

# v29
# gencode_txdb = makeTxDbFromGFF('/s/genomes/human/hg19/gencode29/gencode.v29lift37.annotation.gtf.gz', format='gtf')
# saveDb(gencode_txdb, "/s/genomes/human/hg19/gencode29/gencode.v29lift37.annotation.Db")
gencode_txdb <- loadDb('/s/genomes/human/hg19/gencode29/gencode.v29lift37.annotation.Db')
# Subset to include only canonical chromosomes
std_chr = paste0('chr',c('X','Y','M',1:22))
seqlevels(gencode_txdb) <- std_chr  # Subset the whole txdb object


# Get gene annotation
gtf_or <- rtracklayer::import("/s/genomes/human/hg19/gencode29/gencode.v29lift37.annotation.gtf.gz") %>% as.data.table
gtf_dt <- copy(gtf_or)
gtf_dt <- gtf_dt[type == "gene", .(seqnames, start, end, strand, gene_id, gene_name, gene_type, gene_status)]
gtf_dt <- gtf_dt[seqnames %in% std_chr]
setnames(gtf_dt, "gene_id", "gene_id_full")
gtf_dt <- separate(gtf_dt, "gene_id_full", into = "gene_id", sep = "\\.", remove = F)
head(gtf_dt)

dup_genes <- gtf_dt[duplicated(gtf_dt$gene_name), gene_name] # Get genes that appear at least twice
# Get genes that appear more than twice
repeated_genes <- names(table(gtf_dt[gene_name %in% dup_genes, gene_name])[table(gtf_dt[gene_name %in% dup_genes, gene_name]) > 1])

# rename duplicate gene names
gtf_dt[, N := 1:.N, by = gene_name]
gtf_dt[, gene_name_unique := gene_name]
gtf_dt[N > 1, gene_name_unique := paste(gene_name, N, sep = '_')]
gtf_dt[, N := NULL]

# check if successful
gtf_dt[gene_name %in% repeated_genes]
gtf_dt[duplicated(gtf_dt$gene_id)]  # 41 X-Y paralog genes
dup_genes <- gtf_dt[duplicated(gtf_dt$gene_name), gene_name]
View(gtf_dt[gene_name %in% dup_genes])

fwrite(gtf_dt, "resources/gencode_v29_unique_gene_name.tsv", sep = '\t')


#### NEW ANNOTATION AND EXON EXTRACTION
invert_strand <- function(exons) {
    exons_op <- copy(exons)
    if (class(exons) == "GRanges") {
        strand(exons_op) <- ifelse(strand(exons) == '+', '-', '+')
    } else if (class(exons) %in% c("GRangesList", "CompressedGRangesList")) {
        unlisted <- unlist(exons_op) # create genomic ranges object
        strand(unlisted) <- ifelse(strand(unlisted) == '+', '-', '+')
        exons_op <- relist(unlisted, exons_op) # change back to GRList
    } else {
        message(class(exons))
    }
    exons_op
}

# gtex gencode 19
gtex_op <- invert_strand(exonsBy(gtex_txdb, by = "gene"))
saveRDS(gtex_op, "resources/exons_by_gene_op_v19.rds")
# gencode 29
gencode_op <- invert_strand(exonsBy(gencode_txdb, by = "gene"))
saveRDS(gencode_op, "resources/exons_by_gene_op_v29.rds")

# all genes
# Make sure that the chromosomes names from the bam and annotation files are the same (eg, 1 != chr1)
genes_gr = sort(genes(ucsc_txdb))
genes_en = sort(genes(gtex_txdb))
genes_en = sort(genes(gencode_txdb))


transcripts_en = sort(transcripts(gtex_txdb))
exons_en = sort(exonsBy(gtex_txdb, by = "gene"))

saveRDS(exons_en, "./resources/exons_en.Rds")
exons_en <- readRDS("./resources/exons_en.Rds")


genes_dt = as.data.table(genes_en)
g2 = merge(genes_dt, gene_mapping, by = "gene_id")
saveRDS(g2, "./resources/gencode.v19_with_gene_name.Rds")

exons_en = readRDS("./resources/exons_en.Rds")

genes_gr= with(genes_gr, genes_gr[seqnames %in% std_chr,])
seqlevels(genes_gr) = as.character(seqnames(genes_gr)@values)

genes_gr= with(genes_gr, genes_gr[seqnames %in% std_chr,])
seqlevels(genes_gr) = as.character(seqnames(genes_gr)@values)
# saveRDS(genes_gr, '/s/genomes/human/hg19/ucsc/ucsc_translated_genes_gr.RDS')

## We might need to add _2, _3, ... for repeated gene names
# We can use the following code
genes_dt <- readRDS("./resources/gencode.v19_with_gene_name.Rds")
genes_dt[, N := .N, by = gene_name]
genes_dt[N > 1]
genes_dt[, gene_new_name := gene_name]
genes_dt[N > 1, gene_new_name := paste(gene_name, N, sep = "_")]
genes_dt[, gene_name := gene_new_name]
genes_dt[, N := NULL]
genes_dt[, gene_new_name := NULL]
saveRDS(genes_dt, "./resources/gencode.v19_with_gene_name.Rds")



