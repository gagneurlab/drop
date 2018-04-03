#'---
#' title: Candidate genes from abberrant expression using outrider
#' author: Agne Matuseviciute
#' wb:
#'   input: 
#'     - file_rna_b0: "`sm config['PROC_DATA'] + 'rna_batch0_new_annotation.txt'`"
#'     - file_rna_b1: "`sm config['PROC_DATA'] + 'rna_batch1.txt'`"
#'     - file_rna_b2: "`sm config['PROC_DATA'] + 'rna_batch2_strand_specific.txt'`"
#'     - file_rna_b3: "`sm config['PROC_DATA'] + 'rna_batch3_strand_specific.txt'`"
#'   #output: [
#'     #candidates: "`sm config['PROC_RESULTS'] + 'candidates_outrider.tsv'`"
#'   #]
#' output: 
#'   html_document:
#'     toc_float: yes
#'     toc: yes
#'---
#'

#+ echo=F, warning=F, message=F, results='hide'
source("src/r/config.R")
candidates <- snakemake@output[['candidates']]

file_rna_b0 <- snakemake@input[['file_rna_b0']]
file_rna_b1 <- snakemake@input[['file_rna_b1']]
file_rna_b2 <- snakemake@input[['file_rna_b2']]
file_rna_b3 <- snakemake@input[['file_rna_b3']]
file_disease_info <- file.path(RAWDIR, "gene_info/meta_disease_genes.tsv")

gtf <- "resources/gencode.v19.genes.patched_contigs.gtf.gz" 
gg <- "resources/gencode.v19_with_gene_name.Rds"
sample_anno <- fread("../sample_annotation/Data/sample_annotation.tsv")
omim_dt <- readRDS("/s/project/genetic_diagnosis/resource/omim_dt.Rds")
gn <- readRDS(gg)
ensembl <-gn[,1]
hgnc <-gn[,13]

#' # Read all batches
#+ 
GENE_ANNO <- fread(file_disease_info, na.strings=c('NA',''))
b0<- read.table(file_rna_b0, sep=" ", check.names=F, header=T, row.names=1)
b1<- read.table(file_rna_b1, sep=" ", check.names=F, header=T, row.names=1)
b2<- read.table(file_rna_b2, sep=" ", check.names=F, header=T, row.names=1)
b3<- read.table(file_rna_b3, sep=" ", check.names=F, header=T, row.names=1)

#' # merge only new batches
#+ 
b12 <- merge(b1, b2, by=0, all=FALSE)
b_12 <- b12[,-1]
rownames(b_12) <-  b12[,1]
b123 <- merge(b_12, b3, by=0, all=FALSE)
b_123 <- b123[,-1]
m <- match(b123[,1], hgnc$transcript_name)
rownames(b_123) <- ensembl$gene_id[m]

#' # merge all batches
#+ 
m <- match(rownames(b0), hgnc$transcript_name)
rownames(b0) <- ensembl$gene_id[m]
b0123 <- merge(b0, b_123, by=0, all=FALSE)
b_0123 <- b0123[,-1]
rownames(b_0123) <-  b0123[,1]             

#' # run OUTRIDER, all batches
#+ 
ods_all <- OUTRIDER::OutriderDataSet(countData=b_0123)
ods_all <- OUTRIDER::filterExpression(ods_all, gtfFile=gtf,
         filterGenes=FALSE, savefpkm=TRUE)
OUTRIDER::plotFPKM(ods_all) #set filter genes to FALSE before
ods_all <- OUTRIDER::filterExpression(ods_all, gtfFile=gtf,
         filterGenes=TRUE, savefpkm=TRUE)
ods_all <- OUTRIDER::estimateSizeFactors(ods_all)
ods_all <- OUTRIDER::autoCorrect(ods_all)
#OUTRIDER::plotCountCorHeatmap(ods_all, normalized = TRUE, nCluster = 0, dendrogram="row", names="none", margins=c(1,1,1,1))
ods_all <- OUTRIDER::fit(ods_all)
ods_all <- OUTRIDER::computePvalues(ods_all, alternative="two.sided")
ods_all <- OUTRIDER::computeZscores(ods_all)
#res_all_0.1 <- OUTRIDER::results(ods_all, padj=0.1)
res_all_0.05 <- OUTRIDER::results(ods_all)
OUTRIDER::plotAberrantPerSample(ods_all, padj=0.05)

#' # prepare result tables, all batches. Adding OMIM and sample annotation info
#+ 
m <- match(res_all_0.05$geneID, ensembl$gene_id)
res_all_0.05$hgnc <- hgnc$transcript_name[m]
display_dt <- merge(
    GENE_ANNO[,.(HGNC_GENE_NAME, DISEASE, SUB_CATEGORY,FUNCTION, MIM_NUMBERS)], 
    res_all_0.05, 
    by.x='HGNC_GENE_NAME', 
    by.y='hgnc', 
    all.y=T
)
#omim
disp <- left_join(display_dt, unique(omim_dt[SYMBOL != "", .(SYMBOL, GMIM, PINH, PMIM)]), 
                  by = c("HGNC_GENE_NAME" = "SYMBOL")) %>% as.data.table()
#sample anno
display <- left_join(disp, sample_anno[, .(RNA_ID, KNOWN_MUTATION, CANDIDATE_GENE, COMMENT, SOLVED, 
                               BIOCHEMICAL_DEFECT,  DISEASE, RNA_PERSON)], #CLINICAL_SYMPTOMS,
                  by = c("sampleID" = "RNA_ID")) %>% as.data.table
setnames(display, old = "sampleID", "RNA_ID")
setnames(display, old = "DISEASE.x", "Disease_info")
setnames(display, old = "DISEASE.y", "Disease_sample_anno")


#'  
#'  Full result table, after running OUTRIDER, absolute Z-score > 3, adjusted P-value < 0.05
#' 
#+ echo=F
#head(display_dt)
DT::datatable(
    display, 
    filter='top', 
    rownames = FALSE
)

few_outliers_per_sample <- display[display$AberrantBySample <length(ods_all)*0.001,]

discarded <- display[display$AberrantBySample > length(ods_all)*0.001,c("RNA_ID","RNA_PERSON")]

#'  
#'  Discarded samples with more than 0.1% outliers: 
#' 
#+ echo=F
#head(discarded)
DT::datatable(
    discarded, 
    filter='top', 
    rownames = FALSE
)

#'  
#'  Result table, Samples with less than 0.1% of outliers 
#' 
#+ echo=F
#head(display_dt)
DT::datatable(
    few_outliers_per_sample, 
    filter='top', 
    rownames = FALSE
)

disp_not_solved <-display[display$SOLVED!="TRUE",]

#'  
#'  Samples, which are not solved yet
#' 
#+ echo=F
#head(display_dt)
DT::datatable(
    disp_not_solved, 
    filter='top', 
    rownames = FALSE
)

#' # run OUTRIDER, only new batches
#+ 
ods_123 <- OUTRIDER::OutriderDataSet(countData=b_123)
ods_123 <- OUTRIDER::filterExpression(ods_123, gtfFile=gtf,
         filterGenes=FALSE, savefpkm=TRUE)
OUTRIDER::plotFPKM(ods_123) #set filter genes to FALSE before
ods_123 <- OUTRIDER::filterExpression(ods_123, gtfFile=gtf,
         filterGenes=TRUE, savefpkm=TRUE)
ods_123 <- OUTRIDER::estimateSizeFactors(ods_123)
ods_123 <- OUTRIDER::autoCorrect(ods_123)
#OUTRIDER::plotCountCorHeatmap(ods_123, normalized = TRUE, nCluster = 0, dendrogram="row", names="none", margins=c(1,1))
ods_123 <- OUTRIDER::fit(ods_123)
ods_123 <- OUTRIDER::computePvalues(ods_123, alternative="two.sided")
ods_123 <- OUTRIDER::computeZscores(ods_123)
#res_all_0.1 <- OUTRIDER::results(ods_all, padj=0.1)
res_123_0.05 <- OUTRIDER::results(ods_123)
OUTRIDER::plotAberrantPerSample(ods_123, padj=0.05)


m <- match(res_123_0.05$geneID, ensembl$gene_id)
res_123_0.05$hgnc <- hgnc$transcript_name[m]

display_dt_123 <- merge(
    GENE_ANNO[,.(HGNC_GENE_NAME, DISEASE, SUB_CATEGORY,FUNCTION, MIM_NUMBERS)], 
    res_123_0.05, 
    by.x='HGNC_GENE_NAME', 
    by.y='hgnc', 
    all.y=T
)

#omim
disp_123 <- left_join(display_dt_123, unique(omim_dt[SYMBOL != "", .(SYMBOL, GMIM, PINH, PMIM)]), 
                  by = c("HGNC_GENE_NAME" = "SYMBOL")) %>% as.data.table()
#sample anno
display_123 <- left_join(disp_123, sample_anno[, .(RNA_ID, KNOWN_MUTATION, CANDIDATE_GENE, COMMENT, SOLVED, 
                               BIOCHEMICAL_DEFECT,  DISEASE, RNA_PERSON)], #CLINICAL_SYMPTOMS,
                  by = c("sampleID" = "RNA_ID")) %>% as.data.table
setnames(display_123, old = "sampleID", "RNA_ID")
setnames(display_123, old = "DISEASE.x", "Disease_info")
setnames(display_123, old = "DISEASE.y", "Disease_sample_anno")


few_outliers_per_sample_123 <- display_123[display_123$AberrantBySample <length(ods_all)*0.001,]

discarded_123 <- display_123[display_123$AberrantBySample > length(ods_123)*0.001,c("RNA_ID","RNA_PERSON")]

#'  
#'  Discarded samples with more than 0.1% outliers: 
#' 
#+ echo=F
#head(discarded)
DT::datatable(
    discarded_123, 
    filter='top', 
    rownames = FALSE
)

#'  
#'  Result table, for batches 1, 2, 3, Samples with less than 0.1% of outliers 
#' 
#+ echo=F
#head(display_dt)
DT::datatable(
    few_outliers_per_sample_123, 
    filter='top', 
    rownames = FALSE
)

