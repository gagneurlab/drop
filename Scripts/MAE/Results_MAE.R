#'---
#' title: MAE Results table
#' author: vyepez
#' wb:
#'  input:
#'   - mae_res: '`sm expand(config["PROC_RESULTS"] + "/mae/samples/{id}_res.Rds", id = config["mae_ids"])`'
#'  output:
#'   - res_signif_all: '`sm config["PROC_RESULTS"] + "/mae/MAE_results.Rds"`'
#'   - res_signif_rare: '`sm config["PROC_RESULTS"] + "/mae/MAE_results_rare.Rds"`'
#' output: 
#'   html_document:
#'    code_folding: show
#'    code_download: TRUE
#'---

#+ echo=F
saveRDS(snakemake, 'tmp/mae_res_all.Rds')
# snakemake <- readRDS('tmp/mae_res_all.Rds')
suppressPackageStartupMessages({
    library(data.table)
    library(magrittr)
    library(ggplot2)
    library(cowplot)
    library(tidyr)
})

source("Scripts/_functions/gene_annotation/add_gene_info_cols.R")

#' ## Read all mae files
res <- lapply(snakemake@input$mae_res[51:250], function(m){
    rt <- readRDS(m)
    rt <- rt[padj < .05 & alt_freq > .8]
    return(rt)
}) %>% rbindlist()

res <- separate(res, 'sample', into = c('EXOME_ID', 'RNA_ID'), sep = "-", remove = FALSE)
res[, c('GT', 'as_gt') := NULL]


# Remove mismatches
res <- res[! sample %in% c("EXT_JAP_PT008-103165R", "EXT_JAP_PT875-103207R", "EXT_JAP_PT1146-103229R")]

# Add gene info
setnames(res, "hgncid", "gene_name")
res <- add_all_gene_info(res, dis_genes = F)

# Add sample annotation
sa <- fread(snakemake@config$SAMPLE_ANNOTATION)
res[, gene_name := toupper(gene_name)]
res <- left_join(res, sa[, .(RNA_ID, FIBROBLAST_ID, EXOME_ID, PEDIGREE, KNOWN_MUTATION,
                             CANDIDATE_GENE, BATCH, COMMENT)],
                 by = c("EXOME_ID", "RNA_ID")) %>% as.data.table
setorder(res, padj)

#+ echo=F
res[, aux := paste(chr, pos, REF, ALT, sep = "-")]


# TRUE NAs
res[aux == 'chr1-144951973-T-C', MAX_AF := 0.5]  # misannotated, AF taken from https://gnomad.broadinstitute.org/variant/1-144951973-T-C, PDE4DIP gene
res[aux == 'chr16-87436764-A-G', MAX_AF := 0.79]  # misannotated, AF taken from https://gnomad.broadinstitute.org/variant/16-87436764-A-G, MAP1LC3B gene
res[aux == 'chr1-145112313-T-C', MAX_AF := 0.50]  # misannotated, AF taken from https://gnomad.broadinstitute.org/variant/1-145112313-T-C, SEC22B gene
res[aux == 'chr1-145116833-C-T', MAX_AF := 0.50]  # misannotated, AF taken from https://gnomad.broadinstitute.org/variant/1-145116833-C-T, SEC22B gene
res[aux == 'chr1-145116798-G-A', MAX_AF := 0.50]  # misannotated, AF taken from https://gnomad.broadinstitute.org/variant/1-145116798-G-A, SEC22B gene
res[aux == 'chr1-145116027-T-C', MAX_AF := 0.50]  # misannotated, AF taken from https://gnomad.broadinstitute.org/variant/1-145116027-T-C, SEC22B gene
res[aux == 'chr1-145116679-T-C', MAX_AF := 0.50]  # misannotated, AF taken from https://gnomad.broadinstitute.org/variant/1-145116679-T-C, SEC22B gene
res[aux == 'chr1-145103844-G-A', MAX_AF := 0.50]  # misannotated, AF taken from https://gnomad.broadinstitute.org/variant/1-145103844-G-A, SEC22B gene
res[aux == 'chr1-1581881-T-C', MAX_AF := 0.35]  # misannotated, AF taken from https://gnomad.broadinstitute.org/variant/1-1581881-T-C, CDK11B gene
res[aux == 'chr1-144951959-T-C', MAX_AF := 0.17]  # misannotated, AF taken from https://gnomad.broadinstitute.org/variant/1-144951959-T-C, PDE4DIP gene
res[aux == 'chr1-144854083-G-A', MAX_AF := 0.5]  # misannotated, AF taken from https://gnomad.broadinstitute.org/variant/1-144854083-G-A, PDE4DIP gene
res[aux == 'chr1-144854380-T-C', MAX_AF := 0.5]  # misannotated, AF taken from https://gnomad.broadinstitute.org/variant/1-144854380-T-C, PDE4DIP gene
res[aux == 'chr17-18554918-G-A', MAX_AF := 0.15]  # misannotated, AF taken from https://gnomad.broadinstitute.org/variant/17-18554918-G-A, TBC1D28 gene
res[aux == 'chr17-18554954-A-G', MAX_AF := 0.12]  # misannotated, AF taken from https://gnomad.broadinstitute.org/variant/17-18554954-A-G, TBC1D28 gene
res[aux == 'chr6-29913111-A-G', MAX_AF := 0.17]  # misannotated, AF taken from https://gnomad.broadinstitute.org/variant/6-29913111-A-G, HLA-A gene
res[aux == 'chr6-29913135-G-C', MAX_AF := 0.39]  # misannotated, AF taken from https://gnomad.broadinstitute.org/variant/6-29913135-G-C, HLA-A gene
res[aux == 'chr6-29912754-G-A', MAX_AF := 0.51]  # misannotated, AF taken from https://gnomad.broadinstitute.org/variant/6-29912754-G-A, HLA-A gene
res[aux == 'chr6-32632660-T-G', MAX_AF := 0.19]  # misannotated, AF taken from https://gnomad.broadinstitute.org/variant/6-32632660-T-G, HLA-DQB1 gene
res[aux == 'chr17-36358847-T-A', MAX_AF := 0.75]  # misannotated, AF taken from https://gnomad.broadinstitute.org/variant/17-36358847-T-A, RP11-1407O15.2 gene
res[aux == 'chr10-70392470-G-A', MAX_AF := 0.005]  # misannotated, AF taken from https://gnomad.broadinstitute.org/variant/10-70392470-G-A, TET1 gene
res[aux == 'chr1-21752766-A-G', MAX_AF := 0.48]  # misannotated, AF taken from https://gnomad.broadinstitute.org/variant/1-21752766-A-G, NBPF2P gene
res[aux == 'chr1-21752782-T-G', MAX_AF := 0.34]  # misannotated, AF taken from https://gnomad.broadinstitute.org/variant/1-21752782-T-G, NBPF2P gene
res[aux == 'chr17-20492807-C-G', MAX_AF := 0.03]  # misannotated, AF taken from https://gnomad.broadinstitute.org/variant/17-20492807-C-G, ZSWIM5P2 gene
res[aux == 'chr17-20492893-C-T', MAX_AF := 0.03]  # misannotated, AF taken from https://gnomad.broadinstitute.org/variant/17-20492893-C-T, ZSWIM5P2 gene
res[aux == 'chr17-20492954-C-T', MAX_AF := 0.005]  # misannotated, AF taken from https://gnomad.broadinstitute.org/variant/17-20492954-C-T, ZSWIM5P2 gene
res[aux == 'chr4-189660773-C-A', MAX_AF := 0.12]  # misannotated, AF taken from https://gnomad.broadinstitute.org/variant/4-189660773-C-A, no gene
res[aux == 'chr16-74394671-T-C', MAX_AF := 0.45]  # misannotated, AF taken from https://gnomad.broadinstitute.org/variant/16-74394671-T-C, RP11-252A24.2 gene
res[aux == 'chr17-36361992-C-T', MAX_AF := 0.41]  # misannotated, AF taken from https://gnomad.broadinstitute.org/variant/17-36361992-C-T, RP11-1407O15.2 gene
res[aux == 'chr1-16361916-G-A', MAX_AF := 0.53]  # misannotated, AF taken from https://gnomad.broadinstitute.org/variant/1-16361916-G-A, RP11-5P18.10 gene
res[aux == 'chr1-149290383-G-T', MAX_AF := 0.50]  # misannotated, AF taken from https://gnomad.broadinstitute.org/variant/1-149290383-G-T, RP11-403I13.8 gene
res[aux == 'chr9-41954776-C-G', MAX_AF := 0.67]  # misannotated, AF taken from https://gnomad.broadinstitute.org/variant/9-41954776-C-G, RP11-104G3.2 gene
res[aux == 'chr17-44336943-C-T', MAX_AF := 0.18]  # misannotated, AF taken from https://gnomad.broadinstitute.org/variant/17-44336943-C-T, RP11-104G3.2 gene
res[aux == 'chr2-87355925-G-A', MAX_AF := 0.23]  # misannotated, AF taken from https://gnomad.broadinstitute.org/variant/2-87355925-G-A, AC083899.3 gene
res[aux == 'chr12-118684748-A-C', MAX_AF := 0.72]  # misannotated, AF taken from https://gnomad.broadinstitute.org/variant/12-118684748-A-C, TAOK3 gene
res[aux == 'chr1-16973388-C-T', MAX_AF := 0.61]  # misannotated, AF taken from https://gnomad.broadinstitute.org/variant/1-16973388-C-T, MST1P2 gene
res[aux == 'chr1-16973152-A-G', MAX_AF := 0.67]  # misannotated, AF taken from https://gnomad.broadinstitute.org/variant/1-16973152-A-G, MST1P2 gene
res[aux == 'chr1-145116948-C-T', MAX_AF := 0.50]  # misannotated, AF taken from https://gnomad.broadinstitute.org/variant/1-145116948-C-T, no gene
res[aux == 'chr12-9573224-C-T', MAX_AF := 0.57]  # misannotated, AF taken from https://gnomad.broadinstitute.org/variant/12-9573224-C-T, DDX12P gene
res[aux == 'chr12-31242271-G-A', MAX_AF := 0.81]  # misannotated, AF taken from https://gnomad.broadinstitute.org/variant/12-31242271-G-A, DDX11 gene
res[aux == 'chr1-16976585-C-T', MAX_AF := 0.33]  # misannotated, AF taken from https://gnomad.broadinstitute.org/variant/1-16976585-C-T, MST1P2 gene
res[aux == 'chr1-144854277-A-G', MAX_AF := 0.62]  # misannotated, AF taken from https://gnomad.broadinstitute.org/variant/1-144854277-A-G, PDE4DIP gene
res[aux == 'chr1-144994861-C-T', MAX_AF := 0.6]  # misannotated, AF taken from https://gnomad.broadinstitute.org/variant/1-144994861-C-T, PDE4DIP gene
res[aux == 'chr6-32557504-G-A', MAX_AF := 0.12]  # misannotated, AF taken from https://gnomad.broadinstitute.org/variant/6-32557504-G-A, HLA-DRB1 gene
res[aux == 'ch6-33052743-T-C', MAX_AF := 0.016]  # misannotated, AF taken from https://gnomad.broadinstitute.org/variant/6-33052743-T-C, HLA-DPB1 gene
res[aux == 'chr1-21752830-C-T', MAX_AF := 0.27]  # misannotated, AF taken from https://gnomad.broadinstitute.org/variant/1-21752830-C-T, NBPF2P gene
res[aux == 'chr9-41962740-C-T', MAX_AF := 0.51]  # misannotated, AF taken from https://gnomad.broadinstitute.org/variant/9-41962740-C-T, RP11-204M4.2 gene
res[aux == 'chr7-72711133-C-T', MAX_AF := 0.46]  # misannotated, AF taken from https://gnomad.broadinstitute.org/variant/7-72711133-C-T, POM121B gene
res[aux == 'chr10-46960084-G-A', MAX_AF := 0.46]  # misannotated, AF taken from https://gnomad.broadinstitute.org/variant/10-46960084-G-A, SYT15 gene
res[aux == 'chr17-41466096-A-G', MAX_AF := 0.60]  # misannotated, AF taken from https://gnomad.broadinstitute.org/variant/17-41466096-A-G, LINC00910 gene
res[aux == 'chr17-41466108-G-C', MAX_AF := 0.60]  # misannotated, AF taken from https://gnomad.broadinstitute.org/variant/17-41466108-G-C, LINC00910 gene
res[aux == 'chr1-144676646-T-C', MAX_AF := 0.61]  # misannotated, AF taken from https://gnomad.broadinstitute.org/variant/1-144676646-T-C, no gene
res[aux == 'chr9-67331056-A-T', MAX_AF := 0.54]  # misannotated, AF taken from https://gnomad.broadinstitute.org/variant/9-67331056-A-T, BMS1P10 gene
res[aux == 'chr1-14907-A-G', MAX_AF := 0.51]  # misannotated, AF taken from https://gnomad.broadinstitute.org/variant/1-14907-A-G, WASH7P gene
res[aux == 'chr7-76629626-G-A', MAX_AF := 0.43]  # misannotated, AF taken from https://gnomad.broadinstitute.org/variant/7-76629626-G-A, UPK3B gene
res[aux == 'chr1-16383605-A-G', MAX_AF := 0.63]  # misannotated, AF taken from https://gnomad.broadinstitute.org/variant/1-16383605-A-G, CLCNKB gene
res[aux == 'chr17-21195080-A-C', MAX_AF := 0.59]  # misannotated, AF taken from https://gnomad.broadinstitute.org/variant/17-21195080-A-C, MAP2K3 gene
res[aux == 'chr7-64639752-A-T', MAX_AF := 0.003]  # misannotated, AF taken from https://gnomad.broadinstitute.org/variant/7-64639752-A-T, INTS4L1 gene
res[aux == 'chr22-25574507-C-G', MAX_AF := 0.48]  # misannotated, AF taken from https://gnomad.broadinstitute.org/variant/22-25574507-C-G, KIAA1671 gene
res[aux == 'chr13-25671320-A-G', MAX_AF := 0.12]  # misannotated, AF taken from https://gnomad.broadinstitute.org/variant/13-25671320-A-G, PABPC3 gene, has other lower values
res[aux == 'chr13-25671008-C-T', MAX_AF := 0.16]  # misannotated, AF taken from https://gnomad.broadinstitute.org/variant/13-25671008-C-T, PABPC3 gene, has other lower values
res[aux == 'chr13-25670907-C-A', MAX_AF := 0.15]  # misannotated, AF taken from https://gnomad.broadinstitute.org/variant/13-25670907-C-A, PABPC3 gene, has other lower values



# Variant chr6-57512775-T-G on gene PRIM2 has no annotation on gnomAD and ensembl
# Variant chr6-57513317-T-C on gene PRIM2 has no annotation on gnomAD and ensembl
# Variant chr6-57513248-C-G on gene PRIM2 has no annotation on gnomAD and ensembl
# Variant chr6-57513221-C-T on gene PRIM2 has no annotation on gnomAD and ensembl
# Variant chr6-57513182-A-G on gene PRIM2 has no annotation on gnomAD and ensembl
# Variant chr9-68433537-T-G on gene NA has no annotation on gnomAD and ensembl
# Variant chr2-91809227-A-G on gene NA has no annotation on gnomAD and ensembl
# Variant chr10-47133833-A-G on gene HNRNPA1P33 (non-coding) has no annotation on gnomAD and ensembl
# Variant chr10-47133816-T-A on gene HNRNPA1P33 (non-coding) has no annotation on gnomAD and ensembl
# Variant chrY-14107314-T-C on gene MXRA5P1 (non-coding) has no annotation on gnomAD and ensembl
# Variant chrY-14107068-A-C on gene MXRA5P1 (non-coding) has no annotation on gnomAD and ensembl
# Variant chr7-37960275-C-A on gene EPDR1 has no annotation on gnomAD and ensembl
# Variant chr19-10221620-G-C on gene PPAN has no annotation on gnomAD and ensembl

#'
# Bring gene_name column front
res = cbind(res[, .(gene_name)], res[, -"gene_name"])
res_rare <- res[MAX_AF < .001 | is.na(MAX_AF)]
# tail(sort(table(res_rare[! aux %in%c('chr6-57513221-C-T','chr6-57513248-C-G','chr6-57513317-T-C','chrY-14107314-T-C','chrY-21154426-G-A','chr10-47133833-A-G','chr2-91809227-A-G','chr9-68433537-T-G','chr6-57512775-T-G'),aux])), 10)

saveRDS(res, snakemake@output$res_signif_all)
saveRDS(res_rare, snakemake@output$res_signif_rare)

#' ### Download results tables
write.table(res_rare, "/s/public_webshare/project/genetic_diagnosis/results/MAE_results_rare.tsv", sep = "\t", quote = F, row.names = F)

#' [Download MAE rare results table](https://i12g-gagneurweb.informatik.tu-muenchen.de/project/genetic_diagnosis/results/MAE_results_rare.tsv)
DT::datatable(res_rare, caption = "MAE results", style = 'bootstrap')


#' ## Plots
hist(res$alt_freq, breaks = 20)
hist(res_rare$alt_freq, breaks = 20)

ggplot(res[, .N, by = c('sample', 'BATCH')], aes(BATCH, N)) +
    geom_boxplot()
median(res[, .N, by = sample]$N)

ggplot(res_rare[, .N, by = c('sample', 'BATCH')], aes(BATCH, N)) +
    geom_boxplot()
median(res_rare[, .N, by = sample]$N)

setkey(res, sample, chr, pos, REF, ALT)
setkey(res_rare, sample, chr, pos, REF, ALT)
res[, is_rare := F]
res[res_rare, is_rare := T]

events <- res[, .N, by = c('sample', 'is_rare')]
stat_box_data <- function(y, upper_limit = max(events$N)) {
    data.frame(y = upper_limit, label = paste('median =', round(median(y), 1)))
}
ggplot(events, aes(is_rare, N)) +
    geom_boxplot() +
    stat_summary(fun.data = stat_box_data, geom = "text", vjust = -1) +
    coord_trans(y = 'log10')
