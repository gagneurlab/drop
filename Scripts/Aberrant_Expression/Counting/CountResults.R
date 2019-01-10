#'---
#' title: Analyze Counts
#' author: Michaela Muller
#' wb:
#'  input: 
#'  - counts_v19: '`sm config["PROC_RESULTS"] + "/v19/counts/total_counts.Rds"`'
#'  - counts_v29: '`sm config["PROC_RESULTS"] + "/v29/counts/total_counts.Rds"`'
#'  - filtered_counts: '`sm expand(config["PROC_RESULTS"] + "/{annotation}/outrider/filtered_counts.Rds", annotation=config["ANNOTATIONS"])`'
#' output: 
#'   html_document
#'---

saveRDS(snakemake, "tmp/count_analysis.snakemake")
# snakemake <- readRDS("tmp/count_analysis.snakemake")
suppressPackageStartupMessages({
    library(SummarizedExperiment)
    library(data.table)
})

counts_v19 <- readRDS(snakemake@input$counts_v19)
counts_v29 <- readRDS(snakemake@input$counts_v29)

# ods <- readRDS("/s/project/genetic_diagnosis/processed_results/ods_batches2_3_4_ss.Rds")
gene_counts <- readRDS("/s/project/genetic_diagnosis/processed_data/Rds/batches2_3_4_counts_ss.Rds")
dim(gene_counts)

gencode_v29 <- fread("resources/gencode_v29_unique_gene_name.tsv")
gencode_v19 <- readRDS("resources/gencode.v19_with_gene_name.Rds")


#' List of mito genes
genes_to_check <- c("NDUFAF5", "NDUFAF6",
                    "NDUFA2", "NDUFA7", "NDUFA3", "NDUFA11", "NDUFA13",
                    "DNAJC30", "GCDH", "MOCS1",
                    "MRPL12", "MRPL30", "MRPL38", "MRPL45",
                    "MRPS17", "MRPS21",
                    "MSRB3", "MTG1", "RARS2", "SARS2",
                    "SDHAF2", "SLC25A10", "SLC25A26",
                    "TIMM9", "TIMM10B", "TIMM13", "TIMM23",
                    "TK2", "TOMM5", "TSFM", "ACACA", "ACAD11", "FAHD1", "GATC")

robert_ids <- gencode_v19[gene_name %in% genes_robert, gene_id] %>% unique

#' How many mito genes are in the GTF?


#' How many mito genes have at least 1 read for all samples?

 
#' How many mito genes pass the OUTRIDER filter?

gene_counts_nozero <- gene_counts[rowSums(gene_counts) > 0, ]
dim(gene_counts_nozero)

genes_robert %in% gencodev19$gene_name
genes_robert[!genes_robert %in% gencodev19$gene_name]

robert_ids %in% row.names(gene_counts)
robert_ids %in% row.names(gene_counts_nozero)

robert_ids %in% row.names(ods)

genes_robert %in% row.names(ods)


