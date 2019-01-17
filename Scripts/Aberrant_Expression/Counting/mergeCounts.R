#'---
#' title: Merge the counts for all samples
#' author: Michaela Muller
#' wb:
#'  input: 
#'  - counts: '`sm expand(config["PROC_RESULTS"] + "/{{annotation}}/counts/{sampleID}.Rds", sampleID=config["SAMPLE_IDS"])`'
#'  - gene_annot_dt: "/s/project/genetic_diagnosis/resource/gencode_{annotation}_unique_gene_name.tsv"
#'  output:
#'  - counts_ss: '`sm config["PROC_RESULTS"] + "/{annotation}/counts/total_counts_ss.Rds"`'
#'  - counts_ns: '`sm config["PROC_RESULTS"] + "/{annotation}/counts/total_counts_ns.Rds"`'
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

# Read counts
single_counts <- bplapply(snakemake@input$counts, readRDS)
names(single_counts) <- snakemake@config$SAMPLE_IDS
message(paste("read", length(single_counts), "files"))

# Read annotations
gene_annot_dt <- fread(snakemake@input$gene_annot_dt)
sample_anno <- fread(snakemake@config$SAMPLE_ANNOTATION)
sample_anno <- sample_anno[RNA_ID %in% snakemake@config$SAMPLE_IDS]


merge_counts <- function(counts, gene_annot_dt) {
    total_counts <- do.call(cbind, counts)
    colnames(total_counts) <- gsub(".bam", "", colnames(total_counts))
    # Add gene annotation data (rowData)
    rd <- data.table(gene_id_unique = rownames(total_counts))
    rd <- left_join(rd, gene_annot_dt[,.(gene_id_unique, gene_name_unique, gene_type, gene_status)],
                    by = "gene_id_unique")
    rowData(total_counts) <- rd
    rownames(total_counts) <- rd$gene_name_unique
    total_counts
}

# strand-specific
stranded <- sample_anno[as.logical(IS_RNA_SEQ_STRANDED), RNA_ID]
counts_ss <- merge_counts(single_counts[stranded], gene_annot_dt)
saveRDS(counts_ss, snakemake@output$counts_ss)

# non strand-specific
non_stranded <- sample_anno[!as.logical(IS_RNA_SEQ_STRANDED), RNA_ID]
counts_ns <- merge_counts(single_counts[non_stranded], gene_annot_dt)
saveRDS(counts_ss, snakemake@output$counts_ns)
