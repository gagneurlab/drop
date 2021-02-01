#'---
#' title: Get MAE results
#' author: vyepez, mumichae
#' wb:
#'  log:
#'   - snakemake: '`sm str(tmp_dir / "MAE" / "deseq" / "{vcf}--{rna}.Rds")`'
#'  input:
#'   - mae_counts: '`sm cfg.getProcessedDataDir() + "/mae/allelic_counts/{vcf}--{rna}.csv.gz" `'
#'  output:
#'   - mae_res: '`sm cfg.getProcessedResultsDir() + "/mae/samples/{vcf}--{rna}_res.Rds"`'
#'  type: script
#'---

saveRDS(snakemake, snakemake@log$snakemake)

suppressPackageStartupMessages({
    library(stringr)
    library(tMAE)
})

message("Started with deseq")

# Read mae counts
mae_counts <- fread(snakemake@input$mae_counts, fill=TRUE)
mae_counts[, position := as.numeric(position)]

# Sort by chr
mae_counts[,contig := factor(contig, 
                levels = unique(str_sort(mae_counts$contig, numeric = TRUE)))]
mae_counts <- mae_counts[!grep("Default|opcode", contig)]

print("Running DESeq...")
# Function from MAE pkg
rmae <- DESeq4MAE(mae_counts) ## build test for counting REF and ALT in MAE

### Add AF information from gnomAD
if (snakemake@config$mae$addAF == TRUE) {
    print("Adding gnomAD allele frequencies...")
  
    # obtain the assembly from the config
    genome_assembly <- snakemake@config$genomeAssembly
    rmae <- add_gnomAD_AF(rmae, genome_assembly = genome_assembly,
        max_af_cutoff = snakemake@config$mae$maxAF)
} else {
    rmae[, rare := NA]
}

saveRDS(rmae, snakemake@output$mae_res)
