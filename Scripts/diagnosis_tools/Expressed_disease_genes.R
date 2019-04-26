#'---
#' title: Compare expression of disease genes
#' author: vyepez
#' wb:
#'  input: 
#'   - filtered_v29_ov: '/s/project/genetic_diagnosis/processed_results/v29_overlap/outrider/fib/ods_unfitted.Rds'
#'   - hans_table: "/s/project/mitoMultiOmics/raw_data/gene_info/mitochondrial_disorder_genes_prokisch_mayr_cleaned.tsv"
#'   - script_gene_info: "Scripts/_functions/gene_annotation/add_gene_info_cols.R"
#'   - fresard_genes: "/s/project/genetic_diagnosis/resource/fresard_genes.tsv"
#' output: 
#'   html_document:
#'    code_folding: show
#'    code_download: TRUE
#'---

#+ echo=F
saveRDS(snakemake, "tmp/disgene_comp.snakemake")
# snakemake <- readRDS("tmp/disgene_comp.snakemake")
suppressPackageStartupMessages({
    # library(SummarizedExperiment)
    library(data.table)
    library(tidyr)
    library(ggplot2)
    library(magrittr)
    library(OUTRIDER)
})
source(snakemake@input$script_gene_info)

#' Read tables
ov29 <- readRDS(snakemake@input$filtered_v29_ov)
rd <- rowData(ov29) %>% as.data.table()
rm(ov29)
rd[, gene_id :=  dplyr::first(strsplit(gene_id_unique, split = "\\.")[[1]]), by = 1:nrow(rd)]

hans_table <- fread(snakemake@input$hans_table)
mito_genes <- hans_table[DISEASE =='MITO', HGNC_GENE_NAME]

DIR_ext_genes <- "/s/project/mitoMultiOmics/raw_data/gene_info/external"

hema_genes_dt <- fread(file.path(DIR_ext_genes, "Hematology_genelist_ens.txt"), col.names = 'gene_id')

#+echo=F
####### Clean hema genes #####
hema_genes_dt <- hema_genes_dt[gene_id != "ENSG00000268226"]   # corresponds to ENSG00000130826, DKC1 which is already in the list
hema_genes_dt <- hema_genes_dt[gene_id != "ENSG00000267841"]   # corresponds to ENSG00000102145, GATA1 which is already in the list
hema_genes_dt <- hema_genes_dt[gene_id != "ENSG00000265970"]   # corresponds to ENSG00000010704, HFE which is already in the list
hema_genes_dt <- hema_genes_dt[gene_id != "ENSG00000272852"]   # corresponds to ENSG00000105372, RPS19 which is already in the list
hema_genes_dt <- hema_genes_dt[gene_id != "ENSG00000267912"]   # corresponds to ENSG00000015285, WAS which is already in the list
#######

neuro_genes_dt <- fread(file.path(DIR_ext_genes, "Neurology_genelist.txt"), col.names = 'gene_id')

ophtal_genes_dt <- fread(file.path(DIR_ext_genes, "Ophtalmology_genelist_ens.txt"), col.names = 'gene_id')
##### Clean ophtal genes #######
ophtal_genes_dt <- ophtal_genes_dt[gene_id != "ENSG00000268757"]   # corresponds to ENSG00000101986, ABCD1 which is already in the list
ophtal_genes_dt[gene_id == "ENSG00000166748", gene_id := "ENSG00000273540"]  # corresponds to AGBL1
ophtal_genes_dt <- ophtal_genes_dt[gene_id != "ENSG00000269057"]   # corresponds to ENSG00000102001, CACNA1F which is already in the list
ophtal_genes_dt <- ophtal_genes_dt[gene_id != "ENSG00000262875"]   # corresponds to ENSG00000149260, CAPN5 which is already in the list
ophtal_genes_dt <- ophtal_genes_dt[gene_id != "ENSG00000265594"]   # corresponds to ENSG00000106477, CEP41 which is already in the list
ophtal_genes_dt <- ophtal_genes_dt[! gene_id %in% c("ENSG00000232541","ENSG00000227801", "ENSG00000230930", "ENSG00000206290", "ENSG00000235708", "ENSG00000223699")]   # corresponds to ENSG00000204248, COL11A2 which is already in the list
ophtal_genes_dt <- ophtal_genes_dt[gene_id != "ENSG00000269083"]   # corresponds to ENSG00000131174, COX7B which is already in the list
ophtal_genes_dt <- ophtal_genes_dt[gene_id != "ENSG00000266662"]   # corresponds to ENSG00000186765, FSCN2 which is already in the list
ophtal_genes_dt <- ophtal_genes_dt[gene_id != "ENSG00000264499"]   # corresponds to ENSG00000121634, GJA8 which is already in the list
ophtal_genes_dt[gene_id %in% c("ENSG00000260825", "ENSG00000188888"), gene_id := "ENSG00000277399"]  # corresponds to GPR179
ophtal_genes_dt <- ophtal_genes_dt[gene_id != "ENSG00000269411"]   # corresponds to ENSG00000029993, HMGB3 which is already in the list
ophtal_genes_dt <- ophtal_genes_dt[gene_id != "ENSG00000270408"]   # corresponds to ENSG00000101384, JAG1 which is already in the list
ophtal_genes_dt <- ophtal_genes_dt[gene_id != "ENSG00000263243"]   # corresponds to ENSG00000187242, KRT12 which is already in the list
ophtal_genes_dt <- ophtal_genes_dt[gene_id != "ENSG00000259159"]   # corresponds to ENSG00000235718, MFRP which is already in the list
ophtal_genes_dt <- ophtal_genes_dt[gene_id != "ENSG00000261910"]   # corresponds to ENSG00000137474, MYO7A which is already in the list
ophtal_genes_dt <- ophtal_genes_dt[gene_id != "ENSG00000268281"]   # corresponds to ENSG00000102030, NAA10 which is already in the list
ophtal_genes_dt[gene_id == "ENSG00000031544", gene_id := "ENSG00000278570"]  # corresponds to NR2E3
ophtal_genes_dt <- ophtal_genes_dt[gene_id != "ENSG00000265832"]   # corresponds to ENSG00000131779, PEX11B which is already in the list
ophtal_genes_dt <- ophtal_genes_dt[gene_id != "ENSG00000269666"]   # corresponds to ENSG00000102144, PGK1 which is already in the list
ophtal_genes_dt <- ophtal_genes_dt[gene_id != "ENSG00000265228"]   # corresponds to ENSG00000117360, PRPF3 which is already in the list
ophtal_genes_dt <- ophtal_genes_dt[! gene_id %in% c("ENSG00000273109","ENSG00000271442", "ENSG00000273462", "ENSG00000272932", "ENSG00000273469", 
                                                    "ENSG00000272964", "ENSG00000273376", "ENSG00000273224", "ENSG00000273060")]   # corresponds to ENSG00000105618, PRPF31 which is already in the list
ophtal_genes_dt <- ophtal_genes_dt[gene_id != "ENSG00000107618"]   # corresponds to ENSG00000265203, RBP3 which is already in the list
ophtal_genes_dt <- ophtal_genes_dt[gene_id != "ENSG00000271091"]   # corresponds to ENSG00000102218, RP2 which is already in the list
ophtal_genes_dt <- ophtal_genes_dt[gene_id != "ENSG00000268928"]   # corresponds to ENSG00000100075, SLC25A1 which is already in the list
ophtal_genes_dt <- ophtal_genes_dt[gene_id != "ENSG00000268249"]   # corresponds to ENSG00000126953, TIMM8A which is already in the list
ophtal_genes_dt <- unique(ophtal_genes_dt)
######

rd[, gene_name_unique := toupper(gene_name_unique)]
rd <- add_hans_class(rd, "gene_name_unique", return_all_info = F)
rd[, hema_genes := gene_id %in% hema_genes_dt$gene_id]
rd[, neuro_genes := gene_id %in% neuro_genes_dt$gene_id]
rd[, ophtal_genes := gene_id %in% ophtal_genes_dt$gene_id]

#' Make comparison table
comp_dt <- data.table(category = c("all", "in_annot", "counted", "passedFilter"))
comp_dt[, mito := c(length(mito_genes), rd[MITO_DISEASE_GENE == T, .N], rd[MITO_DISEASE_GENE == T & counted1sample == T, .N], rd[MITO_DISEASE_GENE == T & passedFilter == T, .N])]
comp_dt[, hema := c(nrow(hema_genes_dt), rd[hema_genes == T, .N], rd[hema_genes == T & counted1sample == T, .N], rd[hema_genes == T & passedFilter == T, .N])]
comp_dt[, neuro := c(nrow(neuro_genes_dt), rd[neuro_genes == T, .N], rd[neuro_genes == T & counted1sample == T, .N], rd[neuro_genes == T & passedFilter == T, .N])]
comp_dt[, ophtal := c(nrow(ophtal_genes_dt), rd[ophtal_genes == T, .N], rd[ophtal_genes == T & counted1sample == T, .N], rd[ophtal_genes == T & passedFilter == T, .N])]


mt <- melt(comp_dt, variable.name = 'Disease')
mt[, prop := value / max(value), by = Disease]
mt[, category := factor(category, levels = c("all", "in_annot", "counted", "passedFilter"))]

#' Plot and see intersections
#+ fig.height=9
ggplot(mt, aes(category, prop)) + geom_bar(stat = 'identity', position = 'dodge') + 
    geom_text(aes(label = value),  position = position_dodge(width = .8), vjust = -.1) + 
    theme_bw(base_size = 14) + scale_fill_brewer(palette = "Set2") + facet_wrap(~Disease, ncol = 2)

#' Hema genes not in annot
setdiff(hema_genes_dt$gene_id, rd$gene_id)

#' Ophtal genes not in annot
setdiff(ophtal_genes_dt$gene_id, rd$gene_id)

