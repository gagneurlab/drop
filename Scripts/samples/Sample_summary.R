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

#' ## Read and clean global sample annotation
sa <- sa[!is.na(BATCH)]
sa[, solved := !is.na(KNOWN_MUTATION)]
sa[, sample_type := 'patient']
sa[DISEASE == 'HEALTHY' | FIBROBLAST_ID == 'NHDF', sample_type := 'control']
sa[!is.na(TRANSDUCED_GENE) | TISSUE != 'FIBROBLAST' | GROWTH_MEDIUM != 'GLU', sample_type := 'other']

ggplot(sa, aes(BATCH)) + geom_bar(aes(y = ..count.., fill = sample_type)) + 
    geom_text(aes(label = ..count..), stat = 'count' , vjust = -.5) + 
    theme_bw() + scale_fill_canva(palette="Subdued and proffesional")

ggplot(sa, aes(BATCH)) + geom_bar(aes(y = ..count.., fill = RNA_PERSON)) + 
    theme_bw() + scale_fill_ptol()

#' Total number of samples coming from patients (strand and non strand specific)
sa[sample_type == 'patient', .N]
sa[sample_type == 'patient', table(IS_RNA_SEQ_STRANDED)]   

ggplot(sa[sample_type == 'patient'], aes(BATCH)) + geom_bar(aes(y = ..count.., fill = solved)) + 
    geom_text(aes(label = ..count..), stat = 'count' , vjust = -.5) + 
    theme_bw() + scale_fill_colorblind() + ggtitle('Patients only')

#' ## Download sample annotation with patients only
#+ echo=F
write.table(sa[sample_type == 'patient'][order(BATCH)], "/s/public_webshare/project/genetic_diagnosis/results/sample_anno_rna_patients.txt", sep = "\t", quote = F, row.names = F)
#' [Download rna sample annotation](https://i12g-gagneurweb.informatik.tu-muenchen.de/project/genetic_diagnosis/results/sample_anno_rna_patients.txt)
#'

#' ## Questions to answer
#' ### 1. Sample 62336R appears both in Batch 2 and Batch 4, both are repetition from Laura; which should be included?