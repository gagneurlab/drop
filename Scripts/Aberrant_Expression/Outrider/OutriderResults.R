#'---
#' title: Analyze OUTRIDER object and results
#' author: Michaela Mueller, vyepez
#' wb:
#'  input:
#'   - res_ss: '`sm config["PROC_RESULTS"] + "/v29_overlap/outrider/fib_ss/OUTRIDER_results.tsv"`'
#'   - res_ns: '`sm config["PROC_RESULTS"] + "/v29_overlap/outrider/fib_ns/OUTRIDER_results.tsv"`'
#'  output:
#'   - results: '`sm config["PROC_RESULTS"] + "/v29_overlap/outrider/OUTRIDER_results.tsv"`'
#' output:
#'   html_document:
#'    code_folding: show
#'    code_download: TRUE
#'---

#+ echo=F
saveRDS(snakemake, "tmp/outrider_res_fibs.snakemake")
# snakemake <- readRDS("tmp/outrider_res_fibs.snakemake")

suppressPackageStartupMessages({
    library(ggplot2)
    library(data.table)
    library(dplyr)
    library(cowplot)
    library(OUTRIDER)
})

source("Scripts/_functions/gene_annotation/add_gene_info_cols.R")

#' ## Results
res_ss <- fread(snakemake@input$res_ss)
# res_ss <- fread("/s/project/genetic_diagnosis/processed_results/v29_overlap/outrider/fib_ss/OUTRIDER_results.tsv")
dim(res_ss)
res_ss[, STRANDED := 'Specific']

res_ns <- fread(snakemake@input$res_ns)
# res_ns <- fread("/s/project/genetic_diagnosis/processed_results/v29_overlap/outrider/fib_ns/OUTRIDER_results.tsv")
dim(res_ns)
res_ns[, STRANDED := 'Non Specific']

res <- rbind(res_ss, res_ns)

#' ### How many samples with at least one gene
res[, uniqueN(sampleID), by = STRANDED]


#' ### Download results table
write.table(res, snakemake@output[['results']], quote = F, row.names = F, sep = "\t")
write.table(res, "/s/public_webshare/project/genetic_diagnosis/results/OUTRIDER_results.tsv", sep = "\t", quote = F, row.names = F)

#' [Download OUTRIDER results table](https://i12g-gagneurweb.informatik.tu-muenchen.de/project/genetic_diagnosis/results/OUTRIDER_results.tsv)
DT::datatable(res, caption = "OUTRIDER results", style = 'bootstrap', filter = 'top')

res[tp_sample == T, .N, by = sampleID]
rd <- unique(res[,.(sampleID, AberrantBySample, tp_sample)])
ggplot(rd, aes(tp_sample, AberrantBySample)) + geom_boxplot()

ggplot(rd, aes(tp_sample)) + geom_bar(aes(y = ..count..)) + geom_text(aes(label = ..count..), stat = 'count', vjust = -.5) + 
    labs(x = 'Expression Outlier on causal gene')


res[tp_sample == F,.(sampleID, KNOWN_MUTATION)]

ods_ss <- readRDS("/s/project/genetic_diagnosis/processed_results/v29_overlap/outrider/fib_ss/ods.Rds")
ods_ns <- readRDS("/s/project/genetic_diagnosis/processed_results/v29_overlap/outrider/fib_ns/ods.Rds")
rs <- OUTRIDER::results(ods_ss, all = T)
rs[, abs_z := abs(zScore)]
rn <- OUTRIDER::results(ods_ns, all = T)
rs[padjust < 1, .N]
rs[padjust < .05, .N]
rs[pValue < 1e-4, .N]
rs[pValue < 1e-5, .N]
rn[, abs_z := abs(zScore)]

ru <- unique(res[tp_sample == F,.(sampleID, KNOWN_MUTATION)])
ru[, KNOWN_MUTATION := toupper(KNOWN_MUTATION)]
tl <- lapply(1:nrow(ru), function(i){
    sample = ru$sampleID[i]
    gene = ru$KNOWN_MUTATION[i]
    print(sample)
    if(sample %in% colnames(ods_ss)){
        return(rs[sampleID == sample & geneID == gene])
    } else if(sample %in% colnames(ods_ns)){
        return(rn[sampleID == sample & geneID == gene])
    } else return(NULL)
}) %>% rbindlist()
tl[, FC := round(2^l2fc, 3)]
hist(tl$pValue, breaks = 20)
hist(tl$l2fc)
hist(tl$padjust)
plotVolcano(ods_ns, 'MUC1393')
