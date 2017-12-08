# Rscript
# TODO: Parse protein coding genes from a gene annotation file.
# Make them overlap free for usage with HTSeq-count.
# 
# Author: baderd, vyepez
###############################################################################

suppressPackageStartupMessages(source("src/r/config.R"))

# Contains HGSC to ENCODE gene mapping
gene_mapping = readRDS("./resources/GENCODEv19Mapping.RDS")


# UCSC file
ucsc_txdb= makeTxDbFromGFF('/s/genomes/human/hg19/ucsc/ucsc.translated.gtf', format='gtf')
gencode_txdb= makeTxDbFromGFF('./resources/gencode.v19.genes.patched_contigs.gtf.gz', format='gtf')
gencode_txdb = readRDS('./resources/gencode.v19.genes.patched_contigs.Rds')
 
# only canonical chromosomes
std_chr = paste0('chr',c('X','Y','M',1:22))
  
# all genes
# Make sure that the chromosomes names from the bam and annotation files are the same (eg, 1 != chr1)
genes_gr = sort(genes(ucsc_txdb))
genes_en = sort(genes(gencode_txdb))
transcripts_en = sort(transcripts(gencode_txdb))
exons_en = sort(exons(gencode_txdb))

# Genes annotated in the opposite strand
strand(exons_en[strand(exons_en) == "-",]) <- "*"
strand(exons_en[strand(exons_en) == "+",]) <- "-"
strand(exons_en[strand(exons_en) == "*",]) <- "+"


genes_dt = as.data.table(genes_en)
g2 = merge(genes_dt, gene_mapping, by = "gene_id")
saveRDS(g2, "./resources/gencode.v19_with_gene_name.Rds")
  
seqlevels(exons_en) = paste0("chr", seqnames(exons_en)@values)
seqlevels(exons_en) = gsub("MT", "M", seqnames(exons_en)@values)
exons_en = readRDS("./resources/exons_en.Rds")
  
genes_gr= with(genes_gr, genes_gr[seqnames %in% std_chr,])
seqlevels(genes_gr) = as.character(seqnames(genes_gr)@values)
  
genes_gr= with(genes_gr, genes_gr[seqnames %in% std_chr,])
seqlevels(genes_gr) = as.character(seqnames(genes_gr)@values)
# saveRDS(genes_gr, '/s/genomes/human/hg19/ucsc/ucsc_translated_genes_gr.RDS')
  
# make exon bins for DEXSeq
# Without aggregation, exons that overlap with other exons from different genes are simply skipped
exonic_parts_gr= sort(disjointExons(ucsc_txdb, aggregateGenes=FALSE, includeTranscripts=TRUE))
  
exonic_parts_gr= with(exonic_parts_gr, exonic_parts_gr[seqnames %in% std_chr,])
seqlevels(exonic_parts_gr)= as.character(seqnames(exonic_parts_gr)@values)
#    saveRDS(exonic_parts_gr, '/s/genomes/human/hg19/ucsc/ucsc_translated_disjoint_exons_gr.RDS')


#-----------------------------------------------------------
# BAM files
#-----------------------------------------------------------

# get all prokisch RNA from mitomap
all_rna_ids= unique(SAMPLE_ANNOTATION[RNA_BAM.STAR==T & LAB=='PROKISCH' & TISSUE=='FIBROBLAST', RNA_ID])

# set not counted to all
rna_ids_not_counted = all_rna_ids
rna_ids_not_counted = character()

# get already counted ones
outfile_dir = file.path(DATADIR,'counttable_galignment/rna/')
outfile_prefix = 'ucsc_genes_strict_no_strand_run'
se_already_counted = list.files(outfile_dir, pattern = outfile_prefix, full.names = TRUE)

# substract counted ones if present
if(length(se_already_counted)>0){
  rna_ids_already_counted= unlist( sapply(se_already_counted, function(old_se){
    tmpse = readRDS(old_se)
    colnames(tmpse)
  }
  ), use.names = F
  )
  # update
  rna_ids_not_counted = setdiff(all_rna_ids, rna_ids_already_counted)
}

b_files <- substr(scan("./resources/201708_gusic_rna_seq_samples.tsv", what = character()), 3, 100)
b_files <- substr(scan("./resources/201711_nadel_rna_seq_samples.txt", what = character()), 3, 100)
samples <- vapply(strsplit(b_files, "/"), "[", "", 1)
bf = b_files

# make BAM file collection
rna_bam_list_not_counted= BamFileList(
  sapply(old_samples, get_helmholtz_file, 'rna','BAM.STAR'), # function to get BAM file path
  obeyQname=FALSE, 
  yieldSize= 1e6
)
names(rna_bam_list_not_counted)= rna_ids_not_counted

bamfiles <- BamFileList(file.path(RAWDIR, "helmholtz", bf), yieldSize=2000000)
names(bamfiles) <- samples

# compute size for bam files
bam_size= sort(sapply(bam_files@listData, function(x) file.size(x[['path']])/1e9 ))
summary(bam_size)



#-----------------------------------------------------------
# Define count set
#-----------------------------------------------------------
#'
#' You may want to launch large BAM files with more memory 
#' that is only available on few nodes of your cluster.
#' Therefore reduce the chunk you submit here and 
#' pick i.e. the smallest/largest files only
#' 
chunk_size = 20
idx = names(bam_size)[ 1:min(chunk_size, length(bam_size))]
#idx = names(bam_size)[ I(length(bam_size)-chunk_size):length(bam_size)]

bam_list_to_be_counted = bamfiles[idx]

#-----------------------------------------------------------
# Count parallelized
#-----------------------------------------------------------

# tell SLURM what to do
my_bpparam = register_bplapply_for_clustering(slurm = F, 
                                              workers = 20, 
                                              threads = 2, 
                                              memory = 8000, 
                                              jobname = "count_rna"
)


starttime= Sys.time()
# !!! THIS CALL IS SENSITIVE TO VERSION OF bplapply/biocparallel !!!
library(GenomicAlignments)
se_strand = summarizeOverlaps(
  exons_en, 
  bamfiles, 
  mode = 'IntersectionStrict',
  singleEnd = FALSE,
  ignore.strand = F, # TRUE,
  inter.feature = TRUE, 	# TRUE, reads mapping to multiple features are dropped
  fragments = FALSE, 
  # BPPARAM=NULL) # ,  		# should singletons, reads with unmapped pairs and 
  #   other fragments be included in counting
  BPPARAM = my_bpparam    # SerialParam()
)
message('Processed all fibros in: ', format(Sys.time()- starttime))

# outfile= paste0(outfile_dir, outfile_prefix,'3.rse.RDS')
# saveRDS(se, outfile)

se_strand <- readRDS("/s/project/mitoMultiOmics/processed_expression/counts_batch3.Rds")

# Create exon <-> transcript table
exons_tx = exonsBy(gencode_txdb, by = "tx") %>% unlist %>% as.data.table
setorder(exons_tx, start)
setorder(exons_tx, seqnames)

# Create exon <-> gene table
exons_gene = exonsBy(gencode_txdb, by = "gene") 
exons_gene_dt = exons_gene %>% unlist %>% as.data.table
exons_gene_dt[, gene_name := rep(names(exons_gene), lengths(exons_gene))]
setorder(exons_gene_dt, start)
setorder(exons_gene_dt, seqnames)
setnames(exons_gene_dt, old = "gene_name", new = "gene_id")
exons_gene_dt = merge(exons_gene_dt, gene_mapping, by = "gene_id")
exons_gene_dt <- readRDS("./resources/exons_gene_dt.Rds")

# Seems not to be a difference between tx and gene, so aggregate by gene
identical(transcripts_en$tx_name, genes_en$gene_id)


### Strand specific
se = readRDS(file.path(PROC_DATA, "Rds/se_batch2_strand_specific.Rds"))
se = readRDS(file.path(PROC_DATA, "Rds/se_batch3.Rds"))
head(assay(se))

# Add gene info to exon count
se_ex = assay(se) %>% as.data.table()
se_ex[, exon_id := 1:.N]
se_ex = merge(se_ex, exons_gene_dt, by = "exon_id")


# Aggregate exon data by gene and turn into a matrix
se_gene = se_ex[, lapply(.SD, sum), by = gene_name, .SDcols = colnames(assay(se))]
se_matrix = as.matrix(se_gene[, colnames(assay(se)), with = F])
row.names(se_matrix) = se_gene$gene_name
dim(se_matrix)   # 54,356 rows
se_matrix = se_matrix[rowSums2(se_matrix) > 0, ]
dim(se_matrix)   # ~ 30,000 rows

write.table(se_matrix, file.path(PROC_DATA, "rna_batch3_strand_specific.txt"), quote = F, col.names = T, row.names = T)

### Strand non - specific
se_n = readRDS(file.path(PROC_DATA, "Rds/se_batch1_2_non_strand_specific.Rds"))
head(assay(se_n))

# Add gene info to exon count
se_exn = assay(se_n) %>% as.data.table()
se_exn[, exon_id := 1:.N]
se_exn = merge(se_exn, exons_gene_dt, by = "exon_id")

# Aggregate exon data by gene and turn into a matrix
se_gene_n = se_exn[, lapply(.SD, sum), by = gene_name, .SDcols = colnames(assay(se_n))]
se_matrix_n = as.matrix(se_gene_n[, colnames(assay(se_n)), with = F])
row.names(se_matrix_n) = se_gene_n$gene_name
dim(se_matrix_n)   # 54356 rows
se_matrix_n = se_matrix_n[rowSums2(se_matrix_n) > 0, ]
dim(se_matrix_n)   # 37481 rows

write.table(se_matrix_n, file.path(PROC_DATA, "rna_batch1_2_non_strand_specific.txt"), quote = F, col.names = T, row.names = T)

