#'---
#' title: Proteomics Samples Overview
#' author: vyepez
#' wb:
#'  input: 
#'   - protein_sa: '/s/project/mitoMultiOmics/raw_data/proteome/20180927_robert_proteomics/raw_data/P014_10_template_for_row_wise_normalization_BBM_overview.xlsx'
#'   - protein_mat: '/s/project/mitoMultiOmics/raw_data/proteome/20180927_robert_proteomics/TMT_1_14_row_col_norm_all_wo_dup_wo_contr.xlsx'
#' output: 
#'   html_document:
#'    code_folding: show
#'    code_download: TRUE
#'---

#+ echo=F
saveRDS(snakemake, "tmp/prot_overview.snakemake")
# snakemake <- readRDS("tmp/prot_overview.snakemake")
suppressPackageStartupMessages({
    library(readxl)
    library(data.table)
    library(magrittr)
    library(biomaRt)
    library(ggplot2)
    library(ggthemes)
    library(dplyr)
})

#' ## Read and clean both tables
#' ### Read sample annotation
sa_prot <- read_xlsx(snakemake@input$protein_sa) %>% as.data.table
setnames(sa_prot, "ID sample", "PROTEOME_ID")
sa_prot[, Use := TRUE]
sa_prot[grep("errorenous sample|control sample|technical replicates", Comment), Use := FALSE]

#' The proteomes are classified in different Runs, each of around 8 samples + 2 controls. 3 Runs form a Batch. Batch 4 included replicates only, therefore we won't include it.

sa_prot[Use == TRUE, s_run := 1:.N, by = Run]
sa_prot[, Batch := substr(gsub("mix", "", Run), 1, 1)]
sa_prot[, Batch := paste0("B", Batch)]
sa_prot[Use == FALSE, Batch := NA]
sa_prot[Use == TRUE, s_batch := 1:.N, by = Batch]
sa_prot[, Subgroup := substr(gsub("mix", "", Run), 2, 2)]

#' Valid samples in sample annotation (w/o controls or replicates)
sum(sa_prot$Use)

#' Which samples we don't have in our complete sample annotation
sa <- fread("../sample_annotation/Data/sample_annotation.tsv")
setdiff(sa_prot[Use == T, PROTEOME_ID], sa$PROTEOME_ID)

#' ### Read protein intensities matrix
prot_dt <- read_xlsx(snakemake@input$protein_mat) %>% as.data.table
dim(prot_dt)

# remove outlier samples
prot_dt <- prot_dt[, -c("P33281", "P103213")]
dim(prot_dt)

#' Number of proteins that are 0 in all samples
sum(rowSums(prot_dt[, -"ProteinID"]) == 0)

# remove proteins with all zeroes
prot_dt <- prot_dt[ rowSums(prot_dt[, -"ProteinID"]) > 0]

#' Proteins detected in at least 1 sample
nrow(prot_dt)

# Change to matrix
protein_mat <- as.matrix(prot_dt[, -"ProteinID"])
row.names(protein_mat) <- prot_dt$ProteinID

#' Proteins detected in all samples
sum(apply(protein_mat, 1, function(x) all(x > 0)))

#' ## Intersection plot
op <- do.call(cbind, lapply(unique(sa_prot[Use == T, Batch]), function(b){
    samples_b <- sa_prot[Batch == b & Use == T, PROTEOME_ID]
    samples_b
    apply(protein_mat[, samples_b], 1, function(x) all(x > 0))
}) )
colnames(op) <- unique(sa_prot[Use == T, Batch])

library(UpSetR)
upset_dt <- data.table(Protein_ID = row.names(op), op * 1)
upset(upset_dt, mainbar.y.label = "Proteins detected in all samples of each batch")

#' 
samples_take <- sa_prot[s_run == 1, PROTEOME_ID]
detected_prot <- data.table(PROTEOME_ID = colnames(protein_mat), Proteins_detected = colSums(protein_mat>0))
detected_prot <- left_join(detected_prot, sa_prot[,.(PROTEOME_ID, Batch, Subgroup)]) %>% as.data.table
library(ggbeeswarm)
library(plotly)
#+ 
g <- ggplot(detected_prot, aes(Batch, Proteins_detected, col = Subgroup)) + geom_beeswarm() + 
    theme_bw(base_size = 14) + scale_color_canva()
ggplotly(g)

#' ## Map genes to proteins
# create BiomaRt object for mapping
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

# add HGNC symbol
hgnc_mapping <- data.table(getBM(attributes = c("hgnc_symbol", "uniprotswissprot"),
                                 filters = "uniprotswissprot",
                                 values = row.names(protein_mat),
                                 mart = mart
))

hgnc_mapping <- hgnc_mapping[hgnc_symbol != ""]

#' ### check multiple mappings
# uniprot id wise
hgnc_mapping[duplicated(uniprotswissprot) | duplicated(uniprotswissprot, fromLast=TRUE)]
# hgnc symbol wise
hgnc_mapping[duplicated(hgnc_symbol) | duplicated(hgnc_symbol, fromLast=TRUE)][order(hgnc_symbol)]

dim(protein_mat)
dim(hgnc_mapping)
ods_ss <- readRDS("/s/project/genetic_diagnosis/processed_results/v29_overlap/outrider/ss/ods.Rds")

dim(ods_ss)

library(gplots)
venn(list(transcriptome = row.names(ods_ss), 
          proteome = hgnc_mapping$hgnc_symbol   # already takes unique values
))

#' ## Heatscatter of genes with same proteins
library(LSD)
x <- as.matrix(counts(ods_ss))
heatscatter_genes <- function(gene1, gene2, cor = TRUE){
    heatscatter(x[gene1, ], x[gene2, ], log = 'xy', cor = cor, xlab = gene1, ylab = gene2); abline(0,1)
}
heatpairs_genes <- function(genes, cor = TRUE){
    heatpairs(t(x[genes, ]), log = 'xy', cor = cor)
}

heatscatter_genes("SLX1A", "SLX1B")
heatscatter_genes("RGPD5", "RGPD6")
heatscatter_genes("BCL2L2", "BCL2L2-PABPN1")
heatscatter_genes("SMN1", "SMN2")
heatscatter_genes("MAGED4", "MAGED4B")
heatscatter_genes("AARSD1", "PTGES3L-AARSD1")
# heatscatter_genes("CKMT1A", "CKMT1B")  # all genes not present in ods
# heatscatter_genes("HAB1", "HAB2")    # all genes not present in ods


heatpairs_genes(c("F8A1","F8A2","F8A3"))
heatpairs(t(x[intersect(hgnc_mapping[uniprotswissprot == 'P62805', hgnc_symbol],
                        rownames(ods_ss)), ]), log = 'xy', cor = T)
# heatpairs_genes(c("HIST2H3A","HIST2H3C","HIST2H3D")) ", ]  # all genes not present in ods



