#'---
#' title: Merge the counts for all samples
#' author: Michaela Muller
#' wb:
#'  input: 
#'  - counts: '`sm expand(config["PROC_RESULTS"] + "/counts/overlap/{{annotation}}/{sampleID}_counts.Rds", sampleID=config["SAMPLE_IDS"])`'
#'  - gene_annot_dt: "/s/project/genetic_diagnosis/resource/gencode_{annotation}_unique_gene_name.tsv"
#'  output:
#'  - merged_counts: '`sm config["PROC_RESULTS"] + "/counts/overlap/{annotation}/total_counts.Rds"`'
#'  threads: 30
#'  type: script
#'---

saveRDS(snakemake, "tmp/count_all.snakemake")
# snakemake <- readRDS("tmp/count_all.snakemake")

suppressPackageStartupMessages({
    library(BiocParallel)
    library(SummarizedExperiment)
    library(data.table)
    library(dplyr)
})

register(MulticoreParam(snakemake@threads))

single_counts <- bplapply(snakemake@input$counts, readRDS)
message(paste("read", length(single_counts), "files"))

total_counts <- do.call(cbind, single_counts)
colnames(total_counts) <- gsub(".bam", "", colnames(total_counts))

# Add gene annotation data (rowData)
gene_annot_dt <- fread(snakemake@input$gene_annot_dt)

rd <- data.frame(gene_id_unique = row.names(total_counts))

rd <- left_join(rd, gene_annot_dt[,.(gene_id_unique, gene_name_unique, gene_type, gene_status)],
                                   by = "gene_id_unique")
rowData(total_counts) <- rd

saveRDS(total_counts, snakemake@output$merged_counts)
