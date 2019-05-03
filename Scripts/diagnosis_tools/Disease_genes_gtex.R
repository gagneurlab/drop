#'---
#' title: Disease genes expressed in GTEx
#' author: vyepez
#' wb:
#'  input: 
#'   - gtex_genes: '/s/project/genetic_diagnosis/resource/gtex_expressed_genes.tsv'
#'   - omim_genes: "/s/project/genetic_diagnosis/resource/omim_genes.tsv"
#'   - mito_genes: "/s/project/genetic_diagnosis/resource/mito_disease_genes.tsv"
#'   - fresard_genes: "/s/project/genetic_diagnosis/resource/fresard_genes.tsv"
#'   - neuromuscular_genes: "/s/project/genetic_diagnosis/resource/neuromuscular_genes.tsv"
#' output: 
#'   html_document:
#'    code_folding: show
#'    code_download: TRUE
#'---

#+ echo=F
saveRDS(snakemake, "tmp/disgene_gtex.snakemake")
# snakemake <- readRDS("tmp/disgene_gtex.snakemake")
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(magrittr)
  library(cowplot)
})

#'
#' Read expressed genes in different GTEx tissues
gtex_exp_dt <- fread(snakemake@input$gtex_genes)
gtex_exp_dt[, gene := toupper(gene)]

# Get all genes
all_gtex_genes <- data.table(gene_v19 = toupper(unique(gtex_exp_dt$gene)), DISEASE = 'expressed_genes')
all_gtex_genes[, N := .N]
all_gtex_genes[, ORIGIN := 'GTEx']

# Get all protein coding genes
v19_dt <- fread("/s/project/genetic_diagnosis/resource/gencode_v19_unique_gene_name.tsv")
pc_genes <- v19_dt[gene_type == 'protein_coding', .(gene_v19 = toupper(gene_name_unique))]
pc_genes[, DISEASE := 'protein_coding']
pc_genes[, N := .N]
pc_genes[, ORIGIN := 'v19']

#' Read disease genes tables
omim_dt <- fread(snakemake@input$omim_genes)
mito_genes_dt <- fread(snakemake@input$mito_genes)
fresard_genes_dt <- fread(snakemake@input$fresard_genes)
neuromuscular_genes_dt <- fread(snakemake@input$neuromuscular_genes)

DIS_DT <- rbind(omim_dt, mito_genes_dt, fresard_genes_dt, neuromuscular_genes_dt, 
                all_gtex_genes, pc_genes, use.names = T, fill = T)

#' Combine gtex with disease genes
EXP_DT <- data.table()
for(t in unique(gtex_exp_dt$Tissue_specific)){
  for(d in unique(DIS_DT$DISEASE)){
    dg <- DIS_DT[DISEASE == d, gene_v19]
    prop <- length(intersect(dg, gtex_exp_dt[Tissue_specific == t, gene])) / length(dg)
    dt <- data.table(DISEASE = d, Tissue_specific = t, expressed = prop)
    EXP_DT <- rbind(EXP_DT, dt)
  }
}

#' ## Tissue specific
p <- ggplot(EXP_DT ,aes(DISEASE, Tissue_specific, fill = expressed)) +
  geom_tile() + scale_fill_gradient(low = "white", high = "royalblue4")
p + theme_bw() + theme(axis.text.x = element_text(angle = 90)) 

#' ## Tissue general
gtex_general <- gtex_exp_dt[, .(gene = unique(gene)), by = Tissue_general]
GEN_DT <- data.table()
for(t in unique(gtex_general$Tissue_general)){
  for(d in unique(DIS_DT$DISEASE)){
    dg <- DIS_DT[DISEASE == d, gene_v19]
    prop <- length(intersect(dg, gtex_general[Tissue_general == t, gene])) / length(dg)
    dt <- data.table(DISEASE = d, Tissue_general = t, expressed = prop)
    GEN_DT <- rbind(GEN_DT, dt)
  }
}

p <- ggplot(GEN_DT ,aes(DISEASE, Tissue_general, fill = expressed)) +
  geom_tile() + scale_fill_gradient(low = "white", high = "royalblue4")
p + theme_bw() + theme(axis.text.x = element_text(angle = 90))

gtex_acc <- gtex_general[Tissue_general %in% c("Blood", "Skin", "Muscle")]
setnames(gtex_acc, "Tissue_general", "Tissue")
gtex_acc <- rbind(gtex_acc, data.table(Tissue = "B+S", 
                                       gene = gtex_acc[Tissue %in% c("Blood", "Skin"), gene])) %>% unique

gtex_acc <- rbind(gtex_acc, data.table(Tissue = "B+M", 
                                       gene = gtex_acc[Tissue %in% c("Blood", "Muscle"), gene])) %>% unique

gtex_acc <- rbind(gtex_acc, data.table(Tissue = "B+S+M", 
                                       gene = gtex_acc[Tissue %in% c("Blood", "Skin", "Muscle"), gene])) %>% unique
gtex_acc <- rbind(gtex_acc, data.table(Tissue = "All", 
                                       gene = unique(gtex_general$gene)))

ACC_DT <- data.table()
for(t in unique(gtex_acc$Tissue)){
  for(d in unique(DIS_DT$DISEASE)){
    dg <- DIS_DT[DISEASE == d, gene_v19]
    prop <- length(intersect(dg, gtex_acc[Tissue == t, gene])) / length(dg)
    dt <- data.table(DISEASE = d, Tissue = t, expressed = prop)
    ACC_DT <- rbind(ACC_DT, dt)
  }
}

library(ggthemes)
ACC_DT[, Tissue := factor(Tissue, levels = unique(ACC_DT$Tissue))]
ggplot(ACC_DT, aes(Tissue, expressed, color = DISEASE, group = DISEASE)) + geom_point() + geom_line() + 
  theme_bw() + theme(legend.position = 'top') + 
  scale_color_colorblind()
#' ## TODO
#' 1. Obtain more disease genes
#' 2. Extend to splicing