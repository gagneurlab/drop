#'---
#' title: Create file with file and broken chr
#' author: mumichae
#' wb:
#'  input:
#'   - vcf: '`sm standardFileNames(expand(config["RAW_DATA"] + "/{vcf}/exomicout/paired-endout/stdFilenames/{vcf}.vcf.gz", vcf=vcfs))`'
#'  output:
#'   - debug_list: 'debug_list.tsv'
#'  type: script
#'  threads: 40
#'---

saveRDS(snakemake, 'tmp/debug.Rds')
# snakemake <-  readRDS('tmp/debug.Rds')
suppressPackageStartupMessages({
    library(VariantAnnotation)
    library(data.table)
    library(BiocParallel)
})

register(MulticoreParam(snakemake@threads))
chr_list <- bplapply(snakemake@input$vcf[1:10], function(vcf_file) {
    vcf_obj <- readVcf(vcf_file, "hg19")
    broken_chr <- seqlevels(vcf_obj)[!grepl('chr', seqlevels(vcf_obj))]
    if(length(broken_chr) > 0)
        data.table(vcf_file, broken_chr)
})

fwrite(rbindlist(chr_list), snakemake@output$debug_list)
