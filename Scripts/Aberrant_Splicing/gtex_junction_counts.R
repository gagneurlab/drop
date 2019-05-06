#'---
#' title: 
#' author: vyepez
#' wb:
#'  input:
#'   - gencode_dt: "/s/project/genetic_diagnosis/resource/gencode_v29_unique_gene_name.tsv"
#'   - introns_v29: "/s/project/genetic_diagnosis/processed_results/v29/introns.Rds"
#'  output:
#'   - expressed_junctions_gtex: "/s/project/genetic_diagnosis/resource/gtex_expressed_junctions.tsv"
#'  type: script
#'---

saveRDS(snakemake,  "tmp/counted_junctions.snakemake")
# snakemake <- readRDS("tmp/counted_junctions.snakemake")

suppressPackageStartupMessages({
  library(GenomicFeatures)
  library(GenomicRanges)
  library(data.table)
  library(magrittr)
  library(dplyr)
  devtools::load_all("../FraseR")
})


gencode_dt <- fread(snakemake@input$gencode_dt)
# Read introns / junctions annotation
it <- readRDS(snakemake@input$introns_v29)

DIR_gtex_fraser <- "/s/project/gtex-processed/splicing_map"

gtex_dirs <- list.dirs(file.path(DIR_gtex_fraser, "savedObjects"), full.names = F)
gtex_dirs <- gtex_dirs[grep("-filtered", gtex_dirs)]

# Read ods and obtain expressed genes
DT <- lapply(gtex_dirs, function(gf){
  
  # Read fds and obtain its granges
  fds <- readRDS(file.path(DIR_gtex_fraser, "savedObjects", gf, "fds-object.RDS"))
  
  print(gf)
  gr <- granges(fds, type = 'j')
  gr <- gr[seqnames(gr) %in% c(1:22, "X", "Y", "MT")]
  seqlevels(gr) <- seqlevelsInUse(gr)
  seqlevelsStyle(gr) <- "UCSC"
  
  # Overlap fraser counts with annotated junctions
  ov <- findOverlaps(gr, it, type = 'equal')
  
  # Obtain all counted junctions
  junctions_counted <- it[to(ov)] %>% as.data.table()
  
  max_rep_genes <- max(sapply(junctions_counted$gene_id, length))
  is <- separate(junctions_counted, gene_id, into = paste0("g", 1:max_rep_genes), sep = ",")
  jc <- melt(is, measure.vars = paste0("g", 1:max_rep_genes), value.name = 'gene_id')
  jc <- jc[!is.na(gene_id)]
  jc[, gene_id := gsub('c(', '', gene_id, fixed = T)]
  jc[, gene_id := gsub('\"', '', gene_id, fixed = T)]
  jc[, gene_id := gsub(' ', '', gene_id, fixed = T)]
  jc[, gene_id := gsub(')', '', gene_id, fixed = T)]
  jc[, gene_id := gsub('\n', '', gene_id, fixed = T)]

  jc[, variable := NULL]
  # junctions_counted[, gene_id := as.character(gene_id)]
  jc[, counted_junctions := .N, by = gene_id]
  
  jc <- left_join(jc, gencode_dt[, .(gene_id_unique, gene_name_unique, gene_type, gene_status, N_junctions)], 
                                 by = c('gene_id' = 'gene_id_unique')) %>% as.data.table()
  
  # One row per junction -> one row per gene
  prop_junctions <- unique(jc[, .(gene_id, gene_name_unique, gene_type, gene_status, counted_junctions, N_junctions)])
  setnames(prop_junctions, "N_junctions", "annot_junctions")
  prop_junctions[, prop_junctions := counted_junctions / annot_junctions]
  prop_junctions[, Tissue_specific := gsub("gtex-|-filtered", "", gf)]   # Clean tissue name
  prop_junctions
  
}) %>% rbindlist()


# Read GTEx tissues
gtex_tissues_splicing <- fread("/s/project/genetic_diagnosis/resource/gtex_tissues_splicing.txt")
DT <- left_join(DT, gtex_tissues_splicing, by = 'Tissue_specific') %>% as.data.table()

fwrite(DT, snakemake@output$expressed_junctions_gtex)
