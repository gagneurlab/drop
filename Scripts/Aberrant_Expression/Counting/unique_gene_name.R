#'---
#' title: Create count annotation objects
#' author: mumichae, vyepez
#' wb:
#'  input:
#'   - gtex_gtf: "resources/gencode.v19.genes.patched_contigs.gtf.gz" 
#'   - gencode_gtf: "/s/genomes/human/hg19/gencode29/gencode.v29lift37.annotation.gtf.gz"
#'  output:
#'   - gtex_dt: "/s/project/genetic_diagnosis/resource/gencode_v19_unique_gene_name.tsv"
#'   - gencode_dt: "/s/project/genetic_diagnosis/resource/gencode_v29_unique_gene_name.tsv"
#'  type: script
#'---

saveRDS(snakemake,  "tmp/unique_name.snakemake")
# snakemake <- readRDS("tmp/unique_name.snakemake")
suppressPackageStartupMessages({
    library(BSgenome.Hsapiens.UCSC.hg19)
    library(GenomicFeatures)
    library(GenomicRanges)
    library(data.table)
    library(magrittr)
    library(dplyr)
    library(tidyr)
})

# Function to obtain number of junctions of each gene
obtain_junctions <- function(introns_gr){
  it_dt <- as.data.table(introns_gr)
  max_rep_genes <- max(sapply(it_dt$gene_id, length))
  is <- separate(it_dt, gene_id, into = paste0("g", 1:max_rep_genes), sep = ",")
  jc <- melt(is, measure.vars = paste0("g", 1:max_rep_genes), value.name = 'gene_id')
  jc <- jc[!is.na(gene_id)]
  jc[, gene_id := gsub('c(', '', gene_id, fixed = T)]
  jc[, gene_id := gsub('\"', '', gene_id, fixed = T)]
  jc[, gene_id := gsub(' ', '', gene_id, fixed = T)]
  jc[, gene_id := gsub(')', '', gene_id, fixed = T)]
  jc[, gene_id := gsub('\n', '', gene_id, fixed = T)]

  jc[, variable := NULL]
  
  junctions_dt <- jc[, .(N_junctions = .N), by = gene_id]
  return(junctions_dt)
}


#### v19
# downloaded from https://storage.googleapis.com/gtex_analysis_v7/reference/gencode.v19.genes.v7.patched_contigs.gtf

## Make txdb object
gtex_txdb = makeTxDbFromGFF(snakemake@input$gtex_gtf, format='gtf')
# rename chromosomes from 1 to chr1, ...
seqlevelsStyle(gtex_txdb) <- "UCSC"
message(seqlevels(gtex_txdb))
# saveDb(gtex_txdb, "/s/project/genetic_diagnosis/processed_results/v19/txdb.db")

introns_gtex <- intronicParts(gtex_txdb, linked.to.single.gene.only = FALSE)
saveRDS(introns_gtex, "/s/project/genetic_diagnosis/processed_results/v19/introns.Rds")

## Make dt with gene names
gtex_dt <- rtracklayer::import(snakemake@input$gtex_gtf) %>% as.data.table
gtex_dt <- gtex_dt[type == "transcript", .(seqnames, start, end, width, strand, gene_id, gene_name, gene_type, gene_status)]
gtex_dt[, seqnames := paste0("chr", seqnames)]

## We need to add _2, _3, ... for repeated gene names
gtex_dt[, N := 1:.N, by = gene_name]
gtex_dt[, gene_name_unique := gene_name]
gtex_dt[N > 1, gene_name_unique := paste(gene_name, N, sep = "_")]
gtex_dt[, N := NULL]
gtex_dt[, gene_id_unique := gene_id]

# Add number of junctions
gtex_junctions <- obtain_junctions(introns_gtex)
gtex_dt <- left_join(gtex_dt, gtex_junctions, by = 'gene_id') %>% as.data.table()
gtex_dt[is.na(N_junctions), N_junctions := 0]

fwrite(gtex_dt, snakemake@output$gtex_dt)

#### v29

## Make txdb object
gencode_txdb = makeTxDbFromGFF(snakemake@input$gencode_gtf, format='gtf')
# Subset to include only canonical chromosomes
gencode_txdb <- keepStandardChromosomes(gencode_txdb)
message('v29')
message(seqlevels(gencode_txdb))
# saveDb(gencode_txdb, "/s/project/genetic_diagnosis/processed_results/v29/txdb.db")

introns_gencode <- intronicParts(gencode_txdb, linked.to.single.gene.only = FALSE)
saveRDS(introns_gencode, "/s/project/genetic_diagnosis/processed_results/v29/introns.Rds")

## Make dt with gene names
gtf_dt <- rtracklayer::import(snakemake@input$gencode_gtf) %>% as.data.table
gtf_dt <- gtf_dt[type == "gene", .(seqnames, start, end, strand, gene_id, gene_name, gene_type, gene_status)]
gtf_dt <- gtf_dt[seqnames %in% GenomeInfoDb::standardChromosomes(BSgenome.Hsapiens.UCSC.hg19)]
setnames(gtf_dt, "gene_id", "gene_id_unique")
gtf_dt <- separate(gtf_dt, "gene_id_unique", into = "gene_id", sep = "\\.", remove = F)

dup_genes <- gtf_dt[duplicated(gtf_dt$gene_name), gene_name] # Get genes that appear at least twice
# Get genes that appear more than twice
repeated_genes <- names(table(gtf_dt[gene_name %in% dup_genes, gene_name])[table(gtf_dt[gene_name %in% dup_genes, gene_name]) > 1])

# rename duplicate gene names
gtf_dt[, N := 1:.N, by = gene_name] # warning message
gtf_dt[, gene_name_unique := gene_name]
gtf_dt[N > 1, gene_name_unique := paste(gene_name, N, sep = '_')]
gtf_dt[, N := NULL]

# Add number of junctions
gencode_junctions <- obtain_junctions(introns_gencode)
gtf_dt <- left_join(gtf_dt, gencode_junctions, by = c('gene_id_unique' = 'gene_id')) %>% as.data.table()
gtf_dt[is.na(N_junctions), N_junctions := 0]

fwrite(gtf_dt, snakemake@output$gencode_dt, sep = '\t')


