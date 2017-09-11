# Rscript
# TODO: Parse protein coding genes from a gene annotation file.
# Make them overlap free for usage with HTSeq-count.
# 
# Author: baderd, vyepez
###############################################################################

suppressPackageStartupMessages(source("src/r/config.R"))
REDO_PLOTS=F

#-----------------------------------------------------------
# UCSC file
#-----------------------------------------------------------

if(REDO_PLOTS){
  ucsc_txdb= makeTxDbFromGFF('/s/genomes/human/hg19/ucsc/ucsc.translated.gtf', format='gtf')
  gencode_txdb= makeTxDbFromGFF('/data/ouga/home/ag_gagneur/yepez/workspace/genetic_diagnosis/resources/gencode.v19.genes.patched_contigs.gtf.gz', format='gtf')
  
  # only canonical chromosomes
  std_chr = paste0('chr',c('X','Y','M',1:22))
  
  # all genes
  genes_gr = sort(genes(ucsc_txdb))
  genes_en = sort(genes(gencode_txdb))
  exons_en = sort(exons(gencode_txdb))
  
  genes_gr= with(genes_gr, genes_gr[seqnames %in% std_chr,])
  seqlevels(genes_gr) = as.character(seqnames(genes_gr)@values)
  
  genes_gr= with(genes_gr, genes_gr[seqnames %in% std_chr,])
  seqlevels(genes_gr) = as.character(seqnames(genes_gr)@values)
  #    saveRDS(genes_gr, '/s/genomes/human/hg19/ucsc/ucsc_translated_genes_gr.RDS')
  
  # make exon bins for DEXSeq
  # Without aggregation, exons that overlap with other exons from different genes are simply skipped
  exonic_parts_gr= sort(disjointExons(ucsc_txdb, aggregateGenes=FALSE, includeTranscripts=TRUE))
  
  exonic_parts_gr= with(exonic_parts_gr, exonic_parts_gr[seqnames %in% std_chr,])
  seqlevels(exonic_parts_gr)= as.character(seqnames(exonic_parts_gr)@values)
  #    saveRDS(exonic_parts_gr, '/s/genomes/human/hg19/ucsc/ucsc_translated_disjoint_exons_gr.RDS')
}else{
  genes_gr= readRDS('/s/genomes/human/hg19/ucsc/ucsc_translated_genes_gr.RDS')
  exonic_parts_gr= readRDS('/s/genomes/human/hg19/ucsc/ucsc_translated_disjoint_exons_gr.RDS')
}



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
bf = b_files[1]

# make BAM file collection
rna_bam_list_not_counted= BamFileList(
  sapply(rna_ids_not_counted, get_helmholtz_file, 'rna','BAM.STAR'), # function to get BAM file path
  obeyQname=FALSE, 
  yieldSize= 1e6
)
names(rna_bam_list_not_counted)= rna_ids_not_counted

bamfiles <- BamFileList(file.path(RAWDIR, "helmholtz", bf), yieldSize=2000000)
seqinfo(bamfiles[1])
rna_bam_list_not_counted <- bamfiles

# compute size for bam files
bam_size= sort(sapply(rna_bam_list_not_counted@listData, function(x) file.size(x[['path']])/1e9 ))
summary(bam_size)
#barplot2(bam_size, las=3)



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

bam_list_to_be_counted = rna_bam_list_not_counted[idx]

#-----------------------------------------------------------
# Count parallelized
#-----------------------------------------------------------

# tell SLURM what to do
my_bpparam = register_bplapply_for_clustering(slurm = TRUE, 
                                              workers = length(bam_list_to_be_counted), 
                                              threads = 2, 
                                              memory = 8000, 
                                              jobname = "count_rna_no_strand"
)


starttime= Sys.time()
# !!! THIS CALL IS SENSITIVE TO VERSION OF bplapply/biocparallel !!!
se= summarizeOverlaps(
  genes_gr, 
  bam_list_to_be_counted, 
  mode='IntersectionStrict',
  singleEnd=FALSE,
  ignore.strand= TRUE,
  inter.feature=TRUE, 	# TRUE, reads mapping to multiple features are dropped
  fragments=FALSE, 		# should singletons, reads with unmapped pairs and 
  #   other fragments be included in counting
  BPPARAM = my_bpparam    # SerialParam()
)
message('Processed all fibros in: ', format(Sys.time()- starttime))

outfile= paste0(outfile_dir, outfile_prefix,'3.rse.RDS')
saveRDS(se, outfile)




#-----------------------------------------------------------
# TEST
#-----------------------------------------------------------
if(FALSE){
  sort(log10(sapply(rna_bam_list_not_counted@listData, function(b) file.size(path(b)))))
  
  se= readRDS(outfile)
  se
  summary(assays(se)$counts)
  colData(se)
}

