#'---
#' title: Preprocess Gene Annotations
#' author: mumichae
#' wb:
#'  log:
#'   - snakemake: '`sm str(tmp_dir / "AE" / "{annotation}" / "preprocess.Rds")`'
#'  params:
#'   - assemblyVersion: '`sm cfg.genome.getBSGenomeVersion()`'
#'  input:
#'   - gtf: '`sm lambda wildcards: cfg.genome.getGeneAnnotationFile(wildcards.annotation) `'
#'   - annotate_blacklist: '`sm str(projectDir / ".drop" / "helpers" / "annotate_blacklist.R")`'
#'   - blacklist_19: '`sm str(projectDir / ".drop" / "helpers" / "resource" / "hg19-blacklist.v2.bed.gz")`'
#'   - blacklist_38: '`sm str(projectDir / ".drop" / "helpers" / "resource" / "hg38-blacklist.v2.bed.gz")`'
#'  output:
#'   - txdb: '`sm cfg.getProcessedDataDir() + "/preprocess/{annotation}/txdb.db"`'
#'   - gene_name_mapping: '`sm cfg.getProcessedDataDir() + "/preprocess/{annotation}/gene_name_mapping_{annotation}.tsv"`'
#'   - count_ranges: '`sm cfg.getProcessedDataDir() + 
#'                    "/aberrant_expression/{annotation}/count_ranges.Rds" `'
#'  type: script
#'---

saveRDS(snakemake, snakemake@log$snakemake)

suppressPackageStartupMessages({
  library(GenomicFeatures)
  library(data.table)
  library(rtracklayer)
  library(magrittr)
})

source(snakemake@input$annotate_blacklist)
assemblyVersion <- snakemake@params$assemblyVersion

## Create txdb
txdb <- makeTxDbFromGFF(snakemake@input$gtf)
txdb <- keepStandardChromosomes(txdb)

tmpFile <- tempfile()
saveDb(txdb, tmpFile)
R.utils::copyFile(tmpFile, snakemake@output$txdb, overwrite=TRUE)

# save count ranges
count_ranges <- exonsBy(txdb, by = "gene")
saveRDS(count_ranges, snakemake@output$count_ranges)

# Get a gene annotation table
gtf_dt <- import(snakemake@input$gtf) %>% as.data.table
if (!"gene_name" %in% colnames(gtf_dt)) {
  gtf_dt[gene_name := gene_id]
}
gtf_dt <- gtf_dt[type == 'gene']

if('gene_biotype' %in% colnames(gtf_dt))
  setnames(gtf_dt, 'gene_biotype', 'gene_type')

# Subset to the following columns only
columns <- c('seqnames', 'start', 'end', 'strand', 'gene_id', 'gene_name', 'gene_type', 'gene_status')
columns <- intersect(columns, colnames(gtf_dt))
gtf_dt <- gtf_dt[, columns, with = FALSE]

# Add the blacklist labels for each junction
if(assemblyVersion == 37){
  blacklist_gr <- import(snakemake@input$blacklist_19, format = "BED")
}else{
  blacklist_gr <- import(snakemake@input$blacklist_38, format = "BED")
}

gtf_dt <- addBlacklistLabels(gtf_dt, blacklist_gr, "expression")


# make gene_names unique
gtf_dt[, N := 1:.N, by = gene_name] # warning message
gtf_dt[, gene_name_orig := gene_name]
gtf_dt[N > 1, gene_name := paste(gene_name, N, sep = '_')]
gtf_dt[, N := NULL]



fwrite(gtf_dt, snakemake@output$gene_name_mapping, na = NA)
