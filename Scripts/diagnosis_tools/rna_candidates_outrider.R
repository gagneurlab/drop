#'---
#' title: Candidate genes from abberrant expression using outrider
#' author: Agne Matuseviciute
#' wb:
#'   input: 
#'     - file_rna_b0: "`sm config['PROC_DATA'] + 'rna_batch0_new_annotation.txt'`"
#'     - file_rna_b1: "`sm config['PROC_DATA'] + 'rna_batch1.txt'`"
#'     - file_rna_b2: "`sm config['PROC_DATA'] + 'rna_batch2_strand_specific.txt'`"
#'     - file_rna_b3: "`sm config['PROC_DATA'] + 'rna_batch3_strand_specific.txt'`"
#'   output: [
#'     candidates: "`sm config['PROC_RESULTS'] + 'candidates_outrider.tsv'`"
#'   ]
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
exome_prj <- fread("resources/exome_project_table.csv", na.strings="-")
japan_sampl <- fread("resources/japan_sampl.csv", na.strings="")
remove <- exome_prj$Patient_ID[!is.na(exome_prj$Patient_ID)]
anno_b1 <- fread("resources/MITOMAP_template_b1.csv")
anno_b2 <- fread("resources/MITOMAP_template_b2.csv")
anno_b3 <- fread("resources/batch3.csv")

anno_b2 <-subset(anno_b2, select = -c(IS_RNA_SEQ_STRANDED, V22, V23, V24, V25))
anno_23 <-rbind(anno_b2,anno_b3)

omim_dt <- readRDS("/s/project/genetic_diagnosis/resource/omim_dt.Rds")
mito_carta <-fread("/s/project/mitoMultiOmics/raw_data/gene_info/mito_carta_genes.tsv")
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

b1_names <- colnames(b1)
sample_anno_b0 <- sample_anno[with(sample_anno, !((sample_anno$RNA_ID %in% b1_names))),]
sample_anno_b0 <- sample_anno_b0[,colnames(anno_b1), with=FALSE]
anno_01 <-rbind(sample_anno_b0, anno_b1)
anno_01 <-subset(anno_01, select = -c(IS_RNA_SEQ_STRANDED))

anno_0123 <-rbind(anno_01,anno_23)

############################prepare count data

#' #not strand specific
#+ 
b01 <- merge(b0, b1, by=0, all=FALSE)
b_0_1 <- b01[,-1]
m <- match(b01[,1], hgnc$transcript_name)
rownames(b_0_1) <- ensembl$gene_id[m]

#' #strand specific
#+ 
b23 <- merge(b2, b3, by=0, all=FALSE)
b_2_3 <- b23[,-1]
m <- match(b23[,1], hgnc$transcript_name)
rownames(b_2_3) <- ensembl$gene_id[m]


#' #combine all batches
#+ 
comb <- merge(b_0_1, b_2_3, by=0, all=FALSE)
b_0123 <- comb[,-1]
rownames(b_0123) <- comb[,1]

############################functions

get_outrider_res <- function(count_data, method="BY"){
    ods_all <- OUTRIDER::OutriderDataSet(countData=count_data)
    ods_all <- OUTRIDER::filterExpression(ods_all, gtfFile=gtf,
             filterGenes=TRUE, savefpkm=TRUE)
    ods_all <- OUTRIDER::estimateSizeFactors(ods_all)
    ods_all <- OUTRIDER::autoCorrect(ods_all)
    ods_all <- OUTRIDER::fit(ods_all)
    ods_all <- OUTRIDER::computePvalues(ods_all, alternative="two.sided",  method=method)
    ods_all <- OUTRIDER::computeZscores(ods_all)
    return(ods_all)
}

get_display_table <- function(res, GENE_ANNO, omim_dt, sample_anno){
    #' # prepare result tables, all batches. Adding OMIM and sample annotation info
    #+ 
    m <- match(res$geneID, ensembl$gene_id)
    res$hgnc <- hgnc$transcript_name[m]
    display_dt <- merge(
        GENE_ANNO[,.(HGNC_GENE_NAME, DISEASE, SUB_CATEGORY,FUNCTION, MIM_NUMBERS)], 
        res, 
        by.x='HGNC_GENE_NAME', 
        by.y='hgnc', 
        all.y=T
    )
    #omim
    disp <- left_join(display_dt, unique(omim_dt[SYMBOL != "", .(SYMBOL, GMIM, PINH, PMIM)]), 
                      by = c("HGNC_GENE_NAME" = "SYMBOL")) %>% as.data.table()
    #sample anno
    display <- left_join(disp, sample_anno[, .(RNA_ID, EXOME_ID, GENOME_ID, KNOWN_MUTATION, CANDIDATE_GENE, COMMENT, #SOLVED,
                                   BIOCHEMICAL_DEFECT,  DISEASE, RNA_PERSON, CLINICAL_SYMPTOMS)], #CLINICAL_SYMPTOMS,
                      by = c("sampleID" = "RNA_ID")) %>% as.data.table
    setnames(display, old = "sampleID", "RNA_ID")
    setnames(display, old = "DISEASE.x", "Disease_info")
    setnames(display, old = "DISEASE.y", "Disease_sample_anno")
    return(display)
}


get_final_table <- function(count_data, GENE_ANNO, omim_dt,
                            sample_anno, remove, japan_sampl, mito_carta,
                            b0=NULL, b1=NULL, b2=NULL, b3=NULL){
    ods <-get_outrider_res(count_data, method="BY")
    res_0.05 <- OUTRIDER::results(ods)
    all_res <- OUTRIDER::results(ods, all=T)
    display <-get_display_table(res_0.05, GENE_ANNO, omim_dt, sample_anno)
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
    few_outliers_per_sample <- display[display$AberrantBySample <length(ods)*0.001,]
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
    discarded <- display[display$AberrantBySample > length(ods)*0.001,c("RNA_ID","RNA_PERSON")]
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
    few_outliers_per_sample$KNOWN_MUTATION[few_outliers_per_sample$KNOWN_MUTATION == ""] <- NA
    not_solved <- few_outliers_per_sample[is.na(few_outliers_per_sample$KNOWN_MUTATION)]
    negative_z_all <- not_solved[not_solved$zScore<0,]   
    neutr_all <-all_res[all_res$zScore > -1 & all_res$zScore < 1]
    sample_to_compare <- c()
    sample_to_compare2 <- c()
    sample_to_compare3 <- c()
    for(i in 1:length(negative_z_all$geneID)) {
        sample_to_compare[i] <- neutr_all$sampleID[which(neutr_all$geneID == negative_z_all$geneID[i])[1]]
        sample_to_compare2[i] <- neutr_all$sampleID[which(neutr_all$geneID == negative_z_all$geneID[i])[2]]
        sample_to_compare3[i] <- neutr_all$sampleID[which(neutr_all$geneID == negative_z_all$geneID[i])[3]]
    }
    negative_z_all$compare <- sample_to_compare
    negative_z_all$compare2 <- sample_to_compare2
    negative_z_all$compare3 <- sample_to_compare3
    negative_z_all <- negative_z_all[with(negative_z_all, !((negative_z_all$RNA_ID %in% remove))),]
    if (!is.null(b0)) { 
        negative_z_all$in_b0<- negative_z_all$RNA_ID %in% colnames(b0)
    }
    if (!is.null(b1)) { 
        negative_z_all$in_b1<- negative_z_all$RNA_ID %in% colnames(b1)
    }
    if (!is.null(b2)) { 
        negative_z_all$in_b2<- negative_z_all$RNA_ID %in% colnames(b2)
    }
    if (!is.null(b3)) { 
        negative_z_all$in_b3<- negative_z_all$RNA_ID %in% colnames(b3)
    }
    negative_z_all <- left_join(negative_z_all, japan_sampl[, .(RNA_ID, pedigree, candidate_g_japan)], 
              by = c("RNA_ID" = "RNA_ID")) %>% as.data.table()
    negative_z_all <- inner_join(negative_z_all, mito_carta[, .(HGNC_GENE_NAME, DESCRIPTION, TISSUES, MCARTA_SCORE)], 
                      by = c("HGNC_GENE_NAME" = "HGNC_GENE_NAME")) %>% as.data.table
    negative_z_all[negative_z_all==""]  <- NA
    negative_z_all <- unique(negative_z_all)
    return(negative_z_all)
}


#####################results
    
p_0123 <-get_final_table(b_0123, GENE_ANNO, omim_dt, anno_0123, remove, japan_sampl, mito_carta, b0, b1, b2, b3)    

p_0123_sub <- unique(p_0123_rep[,c("HGNC_GENE_NAME", "geneID","RNA_ID", "EXOME_ID", "GENOME_ID", "PINH",
              "BIOCHEMICAL_DEFECT", "FUNCTION", "zScore", "pValue", "padj_rank",
              "Disease_info","Disease_sample_anno", "PMIM", "GMIM",
              "pedigree", "candidate_g_japan","in_b0","in_b1","in_b2","in_b3"), with=F])


write.table(p_0123, candidates, sep=";", quote=T, row.names=F)

#'  
#'  Final table
#' 
#+ echo=F
#head(display_dt)
DT::datatable(
    p_0123_sub, 
    filter='top', 
    rownames = FALSE
)
