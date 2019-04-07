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
#' Proteomes that we have no transcriptome (including old proteomes)
DT::datatable(sa[!is.na(PROTEOME_ID) & is.na(RNA_ID), .(FIBROBLAST_ID, EXOME_ID, RNA_ID, PROTEOME_ID, PEDIGREE, KNOWN_MUTATION, DISEASE, BATCH)])

#' 61998 and 62343 were mismatches

# Read protein matrix, subset, log-transform and center
protein_gene_mat <- read.csv(snakemake@input$protein_gene_mat) %>% as.matrix
protein_gene_mat[protein_gene_mat < 1e5] <- NA

#' Proteomes that we have no annotation of
setdiff(colnames(protein_gene_mat), sa[,PROTEOME_ID])

rna_prot_dt <- sa[PROTEOME_ID %in% colnames(protein_gene_mat) & TISSUE == 'FIBROBLAST', .(PROTEOME_ID, RNA_ID)]
sa[PROTEOME_ID %in% rna_prot_dt[duplicated(rna_prot_dt$PROTEOME_ID), PROTEOME_ID]]
# Remove one of the replicates
rna_prot_dt = rna_prot_dt[RNA_ID != '110459R']
rna_prot_dt = rna_prot_dt[PROTEOME_ID != 'P62336']  # TODO: check this sample!
rna_prot_dt

protein_gene_mat <- protein_gene_mat[, rna_prot_dt$PROTEOME_ID]
dim(protein_gene_mat)
protein_gene_mat <- log(protein_gene_mat + 1)
protein_gene_mat <- protein_gene_mat - rowMeans2(protein_gene_mat, na.rm = T)
row.names(protein_gene_mat) <- toupper(row.names(protein_gene_mat))


ods <- readRDS("/s/project/genetic_diagnosis/processed_results/v29_overlap/outrider/fib_all/ods.Rds")
ods <- ods[, rna_prot_dt$RNA_ID]
counts <- counts(ods, normalized = F)
counts <- t(t(counts) / sizeFactors(ods))
counts[counts < 50] <- NA
counts <- log(counts + 1) - rowMeans2(log(counts + 1), na.rm = T)

disp <- dispersions(ods)
names(disp) <- rownames(ods)
top_disp_genes <- names(head(sort(disp, decreasing = T), 1000))

common_genes <- intersect(row.names(counts), row.names(protein_gene_mat)) %>% intersect(top_disp_genes)

pm <- protein_gene_mat[common_genes, ]
counts <- counts[common_genes, ]


#' Correlation 
x = sapply(1:ncol(counts), function(i){
    sapply(1:ncol(pm),
    function(j) cor.test(pm[,j], counts[,i], method = 'spearman')$estimate)
    }
 )
hist(x, main = 'Correlation of all RNA - Protein Permutations'); abline(v = 0, col = 'red', lty = 'dashed')
colnames(x) <- rna_prot_dt$RNA_ID
rownames(x) <- rna_prot_dt$PROTEOME_ID
apply(x, 1, which.max)

#' Number of times where the highest correlation was not on the annotated RNA
sum(apply(x, 1, which.max) != 1:nrow(x))

y <- sapply(1:ncol(pm), function(j) cor.test(pm[,j], counts[,j], method = 'spearman')$estimate)
names(y) <- colnames(pm)

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
op <- optim(.3, fr, method = 'Brent', lower = -.5, upper = .8)

plot(mod4, what = "density", data = y, breaks = 15, xlab = "Correlation of RNA - Protein annotated samples")
abline(v = op$par, col = 'red', lty = 'dashed')


#' Example of a mismatch
library(LSD)
#+ fig.width=8, fig.height=8
# heatscatter(counts[, "P104438"], pm[, "P104438"], cor = T, main = 'annotation mismatch'); grid()
# heatscatter(counts[, "P104438"], pm[, "P103200"], cor = T, main = 'actual match'); grid()


#' ## Find possible matches
dt <- as.data.table(melt(x))
setnames(dt, old = c("Var1", "Var2"), c("Proteome_ID", "RNA_ID"))
dt[, max_corr := value == max(value), by = Proteome_ID]
dt[, aux := paste(Proteome_ID, RNA_ID, sep = "-")]
rna_prot_dt[, aux := paste(PROTEOME_ID, RNA_ID, sep = "-")]
dt[, right_annot := aux %in% rna_prot_dt$aux]
dt[, aux := NULL]
DT::datatable(dt[Proteome_ID %in% names(mod$classification[mod$classification == 1]), ][max_corr == T | right_annot == T],  options = list(pageLength = 20))


#' ## CurrFind possible matches


#' ## Plot all correlations of proteomes 
library(plotly)
dt[, Proteome_ID := as.character(Proteome_ID)]
plotlist = list()
for(s in unique(dt$Proteome_ID)){
    g <- ggplot(dt[Proteome_ID == s], aes(reorder(RNA_ID, value), value)) + geom_point(aes(col = right_annot)) + theme_bw() + 
        theme(axis.text.x=element_blank()) + scale_color_calc() + 
        labs(title = s, x = 'Rank', y = 'RNA - Protein Correlation')
    plotlist[[s]] = ggplotly(g)
}

htmltools::tagList(plotlist)
