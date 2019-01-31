#'---
#' title: Samples' summary
#' author: vyepez
#' wb:
#'  input:
#'   - sample_anno: '/s/project/mitoMultiOmics/raw_data/sample_info/SAMPLE_ANNOTATION.tsv'
#'  output:
#' output:
#'   html_document:
#'    code_folding: show
#'    code_download: TRUE
#'---

#+ echo=F
saveRDS(snakemake, "tmp/sample_summary.snakemake")
# snakemake <- readRDS("tmp/sample_summary.snakemake")

suppressPackageStartupMessages({
    library(ggplot2)
    library(ggthemes)
    library(data.table)
})

sa <- fread(snakemake@input$sample_anno)

#'
sa <- sa[!is.na(BATCH)]
sa[, solved := !is.na(KNOWN_MUTATION)]
sa[, sample_type := 'patient']
sa[DISEASE == 'HEALTHY' | FIBROBLAST_ID == 'NHDF', sample_type := 'control']
sa[!is.na(TRANSDUCED_GENE) | TISSUE != 'FIBROBLAST' | GROWTH_MEDIUM != 'GLU', sample_type := 'other']

ggplot(sa, aes(BATCH)) + geom_bar(aes(y = ..count.., fill = sample_type)) + 
    geom_text(aes(label = ..count..), stat = 'count' , vjust = -.5) + 
    theme_bw() + scale_fill_ptol()

ggplot(sa, aes(BATCH)) + geom_bar(aes(y = ..count.., fill = RNA_PERSON)) + 
    theme_bw() + scale_fill_ptol()

ggplot(sa[sample_type == 'patient'], aes(BATCH)) + geom_bar(aes(y = ..count.., fill = solved)) + 
    geom_text(aes(label = ..count..), stat = 'count' , vjust = -.5) + 
    theme_bw() + scale_fill_colorblind() + ggtitle('Patients only')
