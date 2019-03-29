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
  library(OUTRIDER)
})

# Read sample annotation and subset to samples that we have either RNA or Proteome
sa <- fread("../sample_annotation/Data/sample_annotation.tsv")
sa <- sa[! (is.na(RNA_ID) & is.na(PROTEOME_ID))]

#' ## Data Exploration
#' Proteomes that we have no transcriptome
DT::datatable(sa[!is.na(PROTEOME_ID) & is.na(RNA_ID), .(FIBROBLAST_ID, EXOME_ID, RNA_ID, PROTEOME_ID, PEDIGREE, KNOWN_MUTATION, DISEASE, BATCH)])

#' 61998 and 62343 were mismatches

# Read protein matrix, subset, log-transform and center
protein_gene_mat <- read.csv(snakemake@input$protein_gene_mat) %>% as.matrix
protein_gene_mat[protein_gene_mat < 1e5] <- NA
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

ods_ns <- readRDS("/s/project/genetic_diagnosis/processed_results/v29_overlap/outrider/fib_ns/ods.Rds")
ods_ns <- ods_ns[, intersect(colnames(ods_ns), names(prots_ns))]
counts_ns <- counts(ods_ns, normalized = F)
counts_ns <- t(t(counts_ns) / sizeFactors(ods_ns))
counts_ns[counts_ns < 30] <- NA
counts_ns <- log(counts_ns + 1) - rowMeans2(log(counts_ns + 1), na.rm = T)
colnames(counts_ns) <- prots_ns[colnames(counts_ns)]

ods_ss <- readRDS("/s/project/genetic_diagnosis/processed_results/v29_overlap/outrider/fib_ss/ods.Rds")
ods_ss <- ods_ss[, intersect(colnames(ods_ss), names(prots_ss))]
counts_ss <- counts(ods_ss, normalized = F)
counts_ss <- t(t(counts_ss) / sizeFactors(ods_ss))
counts_ss[counts_ss < 50] <- NA
counts_ss <- log(counts_ss + 1) - rowMeans2(log(counts_ss + 1), na.rm = T)
colnames(counts_ss) <- prots_ss[colnames(counts_ss)]

disp <- dispersions(ods_ss)
names(disp) <- rownames(ods_ss)
top_disp_genes <- names(head(sort(disp, decreasing = T), 1000))

common_genes_ns <- intersect(row.names(counts_ns), row.names(pm_ns)) 
common_genes_ss <- intersect(row.names(counts_ss), row.names(pm_ss)) %>% intersect(top_disp_genes)

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

#' Correlation of non strand specific samples
sapply(1:ncol(pm_ns), function(j){
    cor.test(pm_ns[,j], counts_ns[,j], method = 'spearman')$estimate
    })


# register(MulticoreParam(workers = 20))
x = sapply(1:ncol(counts_ss), function(i){
    sapply(1:ncol(pm_ss),
    function(j) cor.test(pm_ss[,j], counts_ss[,i], method = 'spearman')$estimate)
    }
 )
hist(x, main = 'Correlation of all RNA - Protein Permutations'); abline(v = 0, col = 'red', lty = 'dashed')
colnames(x) <- colnames(pm_ss)
rownames(x) <- colnames(pm_ss)
apply(x, 1, which.max)

#' Number of times where the highest correlation was not on the annotated RNA
sum(apply(x, 1, which.max) != 1:nrow(x))

y <- sapply(1:ncol(pm_ss), function(j) cor.test(pm_ss[,j], counts_ss[,j], method = 'spearman')$estimate)
names(y) <- colnames(pm_ss)

#' ## Expectation Maximation to see the classes' separation
library(mclust)
mod <- Mclust(y)
mod4 <- densityMclust(y)
#' Number of groups the algorithm suggests
mod4$G

#' Proteome ids that belong to other cluster
names(mod$classification[mod$classification == 1]) %>% sort

#' Compute the best separation between the classes
mus <- mod4$parameters$mean
sigmas <- sqrt(mod4$parameters$variance$sigmasq)
fr <- function(x) {
    1/sqrt(2*pi*sigmas[1]^2) * exp(-(x-mus[1])^2/(2*sigmas[1]^2)) + 1/sqrt(2*pi*sigmas[2]^2) * exp(-(x-mus[2])^2/(2*sigmas[2]^2))
}
op <- optim(.3, fr, method = 'Brent', lower = -1, upper = 1)

plot(mod4, what = "density", data = y, breaks = 15, xlab = "Correlation of RNA - Protein annotated samples")
abline(v = op$par, col = 'red', lty = 'dashed')


#' Example of a mismatch
library(LSD)
heatscatter(counts_ss[, "P103231"], pm_ss[, "P103231"], cor = T, main = 'annotation mismatch'); grid()
heatscatter(counts_ss[, "P103231"], pm_ss[, "P103230"], cor = T, main = 'actual match'); grid()


#' ## Find possible matches
dt <- as.data.table(melt(x))
setnames(dt, old = c("Var1", "Var2"), c("Proteome_ID", "RNA_ID"))
dt[, max_corr := value == max(value), by = Proteome_ID]
dt[, right_annot := Proteome_ID == RNA_ID]
dt[, match := max_corr == right_annot]
DT::datatable(dt[match == F, ],  options = list(pageLength = 20))


#' ### Plot all correlations of proteomes 
library(plotly)
dt[, Proteome_ID := as.character(Proteome_ID)]
plotlist = list()
for(s in unique(dt$Proteome_ID)){
    g <- ggplot(dt[Proteome_ID == s], aes(reorder(RNA_ID, value), value)) + geom_point(aes(col = right_annot)) + theme_bw() + 
        theme(axis.text.x=element_blank()) + scale_color_calc() + 
        labs(title = s, x = 'Rank', y = 'RNA - Protein Correlation')
    plotlist[[s]] = ggplotly(g)
}

htmltools::tagList(setNames(plotlist, NULL))
