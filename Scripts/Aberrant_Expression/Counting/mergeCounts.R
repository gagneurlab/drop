#'---
#' title: Merge the counts for all samples
#' author: Michaela Muller
#' wb:
#'  input: 
#'    - counts: '`sm lambda wildcards: expand(config["PROC_RESULTS"] + "/{{annotation}}/counts/{sampleID}.Rds", sampleID=config["outrider"][wildcards.dataset])`'
#'    - gene_annot_dt: "/s/project/genetic_diagnosis/resource/gencode_{annotation}_unique_gene_name.tsv"
#'  output:
#'    - counts: '`sm config["PROC_RESULTS"] + "/{annotation}/counts/{dataset}/total_counts.Rds"`'
#'  threads: 30
#'  type: script
#'---
###'  - counts_ss: '`sm expand(config["PROC_RESULTS"] + "/{{annotation}}/counts/{sampleID}.Rds", sampleID=config["outrider"]["fib_ss"])`'
###'  - counts_ns: '`sm expand(config["PROC_RESULTS"] + "/{{annotation}}/counts/{sampleID}.Rds", sampleID=config["outrider"]["fib_ns"])`'
###'  - counts_ss: '`sm config["PROC_RESULTS"] + "/{annotation}/counts/total_counts_ss.Rds"`'
###'  - counts_ns: '`sm config["PROC_RESULTS"] + "/{annotation}/counts/total_counts_ns.Rds"`'

saveRDS(snakemake, "tmp/count_all.snakemake")
# snakemake <- readRDS("tmp/count_all.snakemake")

suppressPackageStartupMessages({
    library(BiocParallel)
    library(SummarizedExperiment)
    library(data.table)
    library(dplyr)
})

register(MulticoreParam(snakemake@threads))

# Read counts
counts <- bplapply(snakemake@input$counts, readRDS)
names(counts) <- snakemake@config$outrider[[snakemake@wildcards$dataset]]
message(paste("read", length(counts), 'files'))

# Read annotations
gene_annot_dt <- fread(snakemake@input$gene_annot_dt)

merge_counts <- function(counts, gene_annot_dt) {
    total_counts <- do.call(cbind, counts)
    colnames(total_counts) <- gsub(".bam", "", colnames(total_counts))
    # Add gene annotation data (rowData)
    rd <- data.table(gene_id_unique = rownames(total_counts))
    rd <- left_join(rd, gene_annot_dt[,.(gene_id_unique, gene_name_unique, gene_type, gene_status)],
                    by = "gene_id_unique")
    rowData(total_counts) <- rd
    total_counts
}

total_counts <- merge_counts(counts, gene_annot_dt)
saveRDS(total_counts, snakemake@output$counts)
