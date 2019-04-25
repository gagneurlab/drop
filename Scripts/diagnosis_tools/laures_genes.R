#'---
#' title: Laure's disease genes
#' author: vyepez
#' wb:
#'  input:
#'   - v19_dt: "/s/project/genetic_diagnosis/resource/gencode_v19_unique_gene_name.tsv"
#'   - v29_dt: "/s/project/genetic_diagnosis/resource/gencode_v29_unique_gene_name.tsv"
#'  output:
#'   - fresard_genes: "/s/project/genetic_diagnosis/resource/fresard_genes.Rds"
#'  type: script
#'---

#+ echo=F
saveRDS(snakemake, "tmp/fresard_genes.snakemake")
# snakemake <- readRDS("tmp/fresard_genes.snakemake")

suppressPackageStartupMessages({
  library(data.table)
  library(tidyr)
  library(magrittr)
  library(dplyr)
})

DIR_ext_genes <- "/s/project/mitoMultiOmics/raw_data/gene_info/external"

# Hematology Genes
hema_genes_dt <- fread(file.path(DIR_ext_genes, "Hematology_genelist_ens.txt"), col.names = 'gene_id')

####### Clean hema genes #####
hema_genes_dt <- hema_genes_dt[gene_id != "ENSG00000268226"]   # corresponds to ENSG00000130826, DKC1 which is already in the list
hema_genes_dt <- hema_genes_dt[gene_id != "ENSG00000267841"]   # corresponds to ENSG00000102145, GATA1 which is already in the list
hema_genes_dt <- hema_genes_dt[gene_id != "ENSG00000265970"]   # corresponds to ENSG00000010704, HFE which is already in the list
hema_genes_dt <- hema_genes_dt[gene_id != "ENSG00000272852"]   # corresponds to ENSG00000105372, RPS19 which is already in the list
hema_genes_dt <- hema_genes_dt[gene_id != "ENSG00000267912"]   # corresponds to ENSG00000015285, WAS which is already in the list
#######

hema_genes_dt[, DISEASE := 'Hematology']

# Neurology genes
neuro_genes_dt <- fread(file.path(DIR_ext_genes, "Neurology_genelist.txt"), col.names = 'gene_id')
neuro_genes_dt[, DISEASE := 'Neurology']

# Ophtalmology genes
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

ophtal_genes_dt[, DISEASE := 'Ophtalmology']

fresard_genes_dt <- rbind(hema_genes_dt, neuro_genes_dt, ophtal_genes_dt)


# Add gene names from v19 and v29
v19_dt <- fread(snakemake@input$v19_dt)
v29_dt <- fread(snakemake@input$v29_dt)
# v19_dt <- fread("/s/project/genetic_diagnosis/resource/gencode_v19_unique_gene_name.tsv")
# v29_dt <- fread("/s/project/genetic_diagnosis/resource/gencode_v29_unique_gene_name.tsv")
v19_dt[, gene_id2 := dplyr::first(unlist(strsplit(gene_id, "\\."))), by = 1:nrow(v19_dt)]

dim(fresard_genes_dt)
fresard_genes_dt <- left_join(fresard_genes_dt, v19_dt[, .(gene_id2, gene_name)], 
                              by = c("gene_id" = "gene_id2")) %>% as.data.table()
dim(fresard_genes_dt)

setnames(fresard_genes_dt, "gene_name", "gene_v19")

fresard_genes_dt <- left_join(fresard_genes_dt, v29_dt[, .(gene_id, gene_name)], by = "gene_id") %>% as.data.table()
dim(fresard_genes_dt)
setnames(fresard_genes_dt, "gene_name", "gene_v29")

fresard_genes_dt[, N := .N, by = DISEASE]
fresard_genes_dt[, ORIGIN := 'Laure Fresard']

fresard_genes_dt[gene_v19 != gene_v29]
fresard_genes_dt[is.na(gene_v19)]

saveRDS(fresard_genes_dt, snakemake@output$fresard_genes)
