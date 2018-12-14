# Author: baderd, vyepez
#########################################################

suppressPackageStartupMessages(source("src/r/config.R"))

# Contains HGSC to ENCODE gene mapping
gene_mapping = readRDS("./resources/GENCODEv19Mapping.RDS")

### UCSC annotation. Already created, just load it.
# ucsc_txdb = makeTxDbFromGFF('/s/genomes/human/hg19/ucsc/ucsc.translated.gtf', format='gtf')
# saveDb(ucsc_txdb, "./resources/ucsc.translated.Db")
ucsc_txdb <- loadDb("./resources/ucsc.translated.Db")

## GTEx annotation. Already created, just load it.
# The gtf file was downloaded from 
# gtex_txdb = makeTxDbFromGFF('./resources/gencode.v19.genes.patched_contigs.gtf.gz', format='gtf')
# gencode_txdb = makeTxDbFromGFF('/s/genomes/human/hg19/gencode29/gencode.v29lift37.annotation.gtf.gz', format='gtf')
gencode_txdb <- loadDb('/s/genomes/human/hg19/gencode29/gencode.v29lift37.annotation.Db')
# saveDb(gencode_txdb, "/s/genomes/human/hg19/gencode29/gencode.v29lift37.annotation.Db")
# saveDb(gtex_txdb, "./resources/gencode.v19.genes.patched_contigs.Db")
gtex_txdb <- loadDb("./resources/gencode.v19.genes.patched_contigs.Db")

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
repeated_genes <- names(table(gtf_dt[gene_name %in% dup_genes, gene_name])[table(gtf_dt[gene_name %in% dup_genes, gene_name]) > 2])

# rename duplicate gene names
for(d in repeated_genes){
    N = nrow(gtf_dt[gene_name == d])
    value = paste(d, 1:N, sep = "_")
    gtf_dt[gene_name == d, gene_name := value]
}
gtf_dt[gene_name %in% repeated_genes]

gtf_dt[duplicated(gtf_dt$gene_id)]  # 41 X-Y paralog genes
dup_genes <- gtf_dt[duplicated(gtf_dt$gene_name), gene_name]
View(gtf_dt[gene_name %in% dup_genes])


#### NEW ANNOTATION AND EXON EXTRACTION
invert_strand <- function(txdb) {
    exons_by_gene <- exonsBy(txdb, by = "gene")
    exons_by_gene_op <- copy(exons_by_gene)
    unlisted <- unlist(exons_by_gene_op)
    strand(unlisted) <- ifelse(strand(unlisted) == '+', '-', '+')
    exons_by_gene_op <- relist(unlisted, exons_by_gene_op)
    exons_by_gene_op
}

# gtex gencode 19
seqlevels(gtex_txdb) <- paste0("chr", seqlevels(gtex_txdb))
seqlevels(gtex_txdb) <- gsub("MT", "M", seqlevels(gtex_txdb))
saveRDS(invert_strand(gtex_txdb), "resources/exons_by_gene_op_v19.rds")
# gencode 29
saveRDS(invert_strand(gencode_txdb), "resources/exons_by_gene_op_v29.rds")

# all genes
# Make sure that the chromosomes names from the bam and annotation files are the same (eg, 1 != chr1)
genes_gr = sort(genes(ucsc_txdb))
genes_en = sort(genes(gtex_txdb))
genes_en = sort(genes(gencode_txdb))


transcripts_en = sort(transcripts(gtex_txdb))
exons_en = sort(exonsBy(gtex_txdb, by = "gene"))

saveRDS(exons_en, "./resources/exons_en.Rds")
exons_en <- readRDS("./resources/exons_en.Rds")

# Genes annotated in the opposite strand
exons_op <- copy(exons_en)
strand(exons_op[strand(exons_op) == "-",]) <- "*"
strand(exons_op[strand(exons_op) == "+",]) <- "-"
strand(exons_op[strand(exons_op) == "*",]) <- "+"
saveRDS(exons_op, "./resources/exons_op.Rds")

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



