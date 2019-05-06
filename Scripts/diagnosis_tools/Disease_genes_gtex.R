#'---
#' title: Disease genes expressed in GTEx
#' author: vyepez
#' wb:
#'  input: 
#'   - gtex_genes: '/s/project/genetic_diagnosis/resource/gtex_expressed_genes.tsv'
#'   - expressed_junctions_gtex: "/s/project/genetic_diagnosis/resource/gtex_expressed_junctions.tsv"
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
  library(ggthemes)
})

#' ## Expression
#' Read expressed genes in different GTEx tissues
gtex_exp_dt <- fread(snakemake@input$gtex_genes)
gtex_exp_dt[, gene := toupper(gene)]

# Get all genes
all_gtex_genes <- data.table(gene_v19 = toupper(unique(gtex_exp_dt$gene)), DISEASE = 'expressed_genes')
all_gtex_genes[, `:=` (N = .N, ORIGIN = 'GTEx')]

# Get all protein coding genes
v19_dt <- fread("/s/project/genetic_diagnosis/resource/gencode_v19_unique_gene_name.tsv")
pc_genes <- v19_dt[gene_type == 'protein_coding', .(gene_v19 = toupper(gene_name_unique))]
pc_genes[, `:=` (DISEASE = 'protein_coding', N = .N, ORIGIN = 'v19')]

#' Read disease genes tables
omim_dt <- fread(snakemake@input$omim_genes)
mito_genes_dt <- fread(snakemake@input$mito_genes)
fresard_genes_dt <- fread(snakemake@input$fresard_genes)
neuromuscular_genes_dt <- fread(snakemake@input$neuromuscular_genes)

DIS_DT <- rbind(omim_dt, mito_genes_dt, fresard_genes_dt, neuromuscular_genes_dt, 
                all_gtex_genes, pc_genes, use.names = T, fill = T)

#' Combine gtex with disease genes
merge_tissues_diseases <- function(EXP_GENES_DT, DISEASE_GENES_DT, tissue_column, gene_version_column = c("gene_v19", "gene_v29")){
  ET <- data.table()
  for(t in unique(EXP_GENES_DT[[tissue_column]])){
    for(d in unique(DISEASE_GENES_DT$DISEASE)){
      dg <- DISEASE_GENES_DT[DISEASE == d, get(gene_version_column)]
      prop <- length(intersect(dg, EXP_GENES_DT[get(tissue_column) == t, gene])) / length(dg)
      dt <- data.table(DISEASE = d, tissue_column = t, expressed = prop)
      ET <- rbind(ET, dt)
    }
  }
  setnames(ET, "tissue_column", tissue_column)
  return(ET)
}

#' ### Tissue specific
EXP_DT <- merge_tissues_diseases(gtex_exp_dt, DIS_DT, "Tissue_specific", "gene_v19")

#+ fig.height=9
p <- ggplot(EXP_DT ,aes(DISEASE, Tissue_specific, fill = expressed)) +
  geom_tile() + scale_fill_gradient(low = "white", high = "royalblue4")
p + theme_bw() + theme(axis.text.x = element_text(angle = 90)) 

#' ### Tissue general
gtex_general <- gtex_exp_dt[, .(gene = unique(gene)), by = Tissue_general]
GEN_DT <- merge_tissues_diseases(gtex_general, DIS_DT, "Tissue_general", "gene_v19")

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

ACC_DT <- merge_tissues_diseases(gtex_acc, DIS_DT, "Tissue", "gene_v19")

ACC_DT[, Tissue := factor(Tissue, levels = unique(ACC_DT$Tissue))]
ggplot(ACC_DT, aes(Tissue, expressed, color = DISEASE, group = DISEASE)) + geom_point() + geom_line() + 
  theme_bw() + theme(legend.position = 'top') + 
  scale_color_colorblind()

#' ## Splicing
#' Read expressed junctions genes in different GTEx tissues
gtex_junc_dt <- fread(snakemake@input$expressed_junctions_gtex)
gtex_junc_dt <- fread("/s/project/genetic_diagnosis/resource/gtex_expressed_junctions.tsv")
setnames(gtex_junc_dt, "gene_name_unique", "gene")
gtex_junc_dt[, gene := toupper(gene)]

# Get all genes
all_gtex_genes <- data.table(gene_v29 = unique(gtex_junc_dt$gene), DISEASE = 'all_genes_cov_junctions')
all_gtex_genes[, `:=` (N = .N, ORIGIN = 'GTEx')]

# Get all protein coding genes
v29_dt <- fread("/s/project/genetic_diagnosis/resource/gencode_v29_unique_gene_name.tsv")
pc_genes <- v29_dt[gene_type == 'protein_coding', .(gene_v29 = toupper(gene_name_unique))]
pc_genes[, `:=` (DISEASE = 'protein_coding', N = .N, ORIGIN = 'v29')]

DIS_JUN_DT <- rbind(omim_dt, mito_genes_dt, fresard_genes_dt, neuromuscular_genes_dt, 
                all_gtex_genes, pc_genes, use.names = T, fill = T)

#' ### Tissue specific
JUNC_DT <- merge_tissues_diseases(gtex_junc_dt, DIS_JUN_DT, "Tissue_specific", "gene_v29")

#+ fig.height=9
p <- ggplot(JUNC_DT ,aes(DISEASE, Tissue_specific, fill = expressed)) +
  geom_tile() + scale_fill_gradient(low = "white", high = "royalblue4")
p + theme_bw() + theme(axis.text.x = element_text(angle = 90)) 

#' ### Tissue general
gtex_junc_general <- gtex_junc_dt[, .(counted_junctions = max(counted_junctions),
                                      annot_junctions = unique(annot_junctions),
                                      prop_junctions = max(prop_junctions)), by = .(gene, Tissue_general)]

GEN_JUNC_DT <- merge_tissues_diseases(gtex_junc_general, DIS_JUN_DT, "Tissue_general", "gene_v29")

p <- ggplot(GEN_JUNC_DT ,aes(DISEASE, Tissue_general, fill = expressed)) +
  geom_tile() + scale_fill_gradient(low = "white", high = "royalblue4")
p + theme_bw() + theme(axis.text.x = element_text(angle = 90)) 

#' ### % of expressed junctions
ggplot(gtex_junc_dt, aes(reorder(Tissue_specific, prop_junctions, FUN = median), prop_junctions)) + geom_boxplot() + theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) + labs(x = 'Tissue')

gj <- gtex_junc_general[, .(m = median(prop_junctions)), by = Tissue_general]
ggplot(gtex_junc_general, aes(reorder(Tissue_general, prop_junctions, FUN = median), 
                              prop_junctions)) + geom_boxplot() + theme_bw() + 
  theme(axis.text.x = element_text(angle = 90)) + labs(x = 'Tissue')

