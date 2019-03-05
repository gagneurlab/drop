#'---
#' title: Protein - RNA correlation
#' author: vyepez
#' wb:
#'  input: 
#'   - protein_gene_mat: '/s/project/mitoMultiOmics/raw_data/proteome/protein_mat_gene_names.txt'
#' output: 
#'   html_document:
#'    code_folding: show
#'    code_download: TRUE
#'---

#+ echo=F
saveRDS(snakemake, "tmp/prot_rna.snakemake")
# snakemake <- readRDS("tmp/prot_rna.snakemake")
suppressPackageStartupMessages({
    library(data.table)
    library(magrittr)
    library(ggplot2)
    library(ggthemes)
    library(dplyr)
    library(matrixStats)
})


sa <- fread("../sample_annotation/Data/sample_annotation.tsv")
sa <- sa[! (is.na(RNA_ID) & is.na(PROTEOME_ID))]
sa[, IS_RNA_SEQ_STRANDED := NULL]
sa[BATCH %in% c("Batch_0", "Batch_1"), IS_RNA_SEQ_STRANDED := F]
sa[BATCH %in% c("Batch_2", "Batch_3", "Batch_4", "Batch_5", "Batch_6"), IS_RNA_SEQ_STRANDED := T]

#'
#' Proteomes that we have no transcriptome
DT::datatable(sa[!is.na(PROTEOME_ID) & is.na(RNA_ID), .(FIBROBLAST_ID, EXOME_ID, RNA_ID, PROTEOME_ID, PEDIGREE, KNOWN_MUTATION, DISEASE, BATCH)])

#' 61998 and 62343 were mismatches

protein_gene_mat <- read.csv(snakemake@input$protein_gene_mat) %>% as.matrix
# protein_gene_mat <- read.csv('/s/project/mitoMultiOmics/raw_data/proteome/protein_mat_gene_names.txt') %>% as.matrix
protein_gene_mat[protein_gene_mat < 10000] <- NA
protein_gene_mat <- log(protein_gene_mat + 1)
protein_gene_mat <- protein_gene_mat - rowMeans2(protein_gene_mat, na.rm = T)
row.names(protein_gene_mat) <- toupper(row.names(protein_gene_mat))

#' Proteomes that we have no annotation of
setdiff(colnames(protein_gene_mat), sa[,PROTEOME_ID])

prots_ns <- sa[!is.na(PROTEOME_ID) & IS_RNA_SEQ_STRANDED == F, PROTEOME_ID]
names(prots_ns) <- sa[!is.na(PROTEOME_ID) & IS_RNA_SEQ_STRANDED == F, RNA_ID]

prots_ss <- sa[!is.na(PROTEOME_ID) & IS_RNA_SEQ_STRANDED == T, PROTEOME_ID]
names(prots_ss) <- sa[!is.na(PROTEOME_ID) & IS_RNA_SEQ_STRANDED == T, RNA_ID]


pm_ns <- protein_gene_mat[, intersect(colnames(protein_gene_mat), prots_ns)]
dim(pm_ns)
pm_ss <- protein_gene_mat[, intersect(colnames(protein_gene_mat), prots_ss)]
dim(pm_ss)

ods_ns <- readRDS("/s/project/genetic_diagnosis/processed_results/v29_overlap/outrider/ns/ods.Rds")
ods_ns <- ods_ns[, intersect(colnames(ods_ns), names(prots_ns))]
counts_ns <- counts(ods_ns, normalized = F)
counts_ns <- t(t(counts_ns) / sizeFactors(ods_ns))
counts_ns[counts_ns < 30] <- NA
counts_ns <- log(counts_ns + 1) - rowMeans2(log(counts_ns + 1), na.rm = T)
colnames(counts_ns) <- prots_ns[colnames(counts_ns)]

ods_ss <- readRDS("/s/project/genetic_diagnosis/processed_results/v29_overlap/outrider/ss/ods.Rds")
ods_ss <- ods_ss[, intersect(colnames(ods_ss), names(prots_ss))]
counts_ss <- counts(ods_ss, normalized = F)
counts_ss <- t(t(counts_ss) / sizeFactors(ods_ss))
counts_ss[counts_ss < 50] <- NA
counts_ss <- log(counts_ss + 1) - rowMeans2(log(counts_ss + 1), na.rm = T)
colnames(counts_ss) <- prots_ss[colnames(counts_ss)]

common_genes_ns <- intersect(row.names(counts_ns), row.names(pm_ns))
common_genes_ss <- intersect(row.names(counts_ss), row.names(pm_ss))

pm_ns <- pm_ns[common_genes_ns, ]
counts_ns <- counts_ns[common_genes_ns, ]

pm_ss <- pm_ss[common_genes_ss, ]
counts_ss <- counts_ss[common_genes_ss, ]


common_ids_ns <- intersect(colnames(counts_ns), colnames(pm_ns))
common_ids_ss <- intersect(colnames(counts_ss), colnames(pm_ss))

pm_ns <- pm_ns[, common_ids_ns]
counts_ns <- counts_ns[, common_ids_ns]

pm_ss <- pm_ss[, common_ids_ss]
counts_ss <- counts_ss[, common_ids_ss]


sapply(1:ncol(pm_ns), function(j){
    cor.test(pm_ns[,j], counts_ns[,j], method = 'spearman')$estimate
    })


register(MulticoreParam(workers = 20))
x = bplapply(1:ncol(counts_ss), function(i){
    sapply(1:ncol(pm_ss),
    function(j) cor.test(pm_ss[,j], counts_ss[,i], method = 'spearman')$estimate)
    }
 )
hist(x, main = 'Correlation of all RNA - Protein Permutations')
colnames(x) <- colnames(pm_ss)
rownames(x) <- colnames(pm_ss)
apply(x, 1, which.max)

#' Number of times where the highest correlation was not on the annotated RNA
sum(apply(x, 1, which.max) != 1:nrow(x))

y <- sapply(1:ncol(pm_ss), function(j) cor.test(pm_ss[,j], counts_ss[,j], method = 'spearman')$estimate)
names(y) <- colnames(pm_ss)
y

boxplot(as.numeric(x),y)

library(LSD)
heatscatter(counts_ss[, "P103231"], pm_ss[, "P103231"], cor = T)

cor.test(pm_ss[, "P103231"], counts_ss[, "P103231"], method = 'spearman')$estimate

dt <- as.data.table(melt(x))
dt[, m1 := value == max(value), by = Var1]
dt[, m2 := Var1 == Var2]
dt[, match := m1 == m2]
dt[match == F, ]

library(plotly)
dt[, Var1 := as.character(Var1)]
plotlist = list()
for(s in unique(dt$Var1)){
    g <- ggplot(dt[Var1 == s], aes(reorder(Var2, value), value)) + geom_point(aes(col = m2)) + theme_bw() + 
        theme(axis.text.x=element_blank()) + scale_color_calc() + 
        labs(title = s, x = 'Rank', y = 'RNA - Protein Correlation')
    plotlist[[s]] = ggplotly(g)
}

htmltools::tagList(setNames(plotlist, NULL))

