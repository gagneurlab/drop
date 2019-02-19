#'---
#' title: MAE cascade plot
#' author: vyepez, mumichae
#' wb:
#'  input:
#'   - out_mae: 'Output/mae.done'
#'  output:
#'  threads: 40
#' output: 
#'   html_document:
#'    code_folding: show
#'    code_download: TRUE
#'---

#+ echo=F
saveRDS(snakemake, 'tmp/cascade_mae.Rds')
# snakemake <- readRDS('tmp/cascade_mae.Rds')
suppressPackageStartupMessages({
    devtools::load_all("../mae/")
    library(ggplot2)
    library(plotly)
    library(magrittr)
})

#' ## Read all mae files
maes <- list.files("/s/project/genetic_diagnosis/processed_data/mae/")
min_cov <- 10

BPPARAM = MulticoreParam(snakemake@threads, progressbar=TRUE)

#+ results = F
dt_mae <- bplapply(maes, function(m){
    mae_table <- readRDS(file.path("/s/project/genetic_diagnosis/processed_data/mae/", m))
    ntotal <- length(mae_table)
    
    # subset for heterozygous mutation only
    mae_table$GT <- as.character(mae_table$GT)
    hetGT <- grepl('0[|/]1|1[|/]0', mae_table$GT, perl=TRUE)
    mae_table <- mae_table[hetGT]
    nhet <- length(mae_table)

    # subset for enough coverage
    mae_table$GQ <- as.integer(mae_table$GQ)
    goodCov <- mae_table$coverage >= min_cov
    mae_table <- mae_table[goodCov]
    nreads <- length(mae_table)
    
    rm(mae_table)
    return(c(ntotal, nhet, nreads))
    }, BPPARAM = BPPARAM
)

#'
names(dt_mae) <- maes

#' ## Get a data.table
m_dt <- melt(data.table(type = c("all", "het", "enough_counts"), as.data.table(dt_mae)), variable.name = 'sample')
m_dt[, sample := as.character(sample)]
m_dt[, sample := head(unlist(strsplit(sample, split = "\\.")), 1), by = 1:nrow(m_dt)]
m_dt[, type := factor(type, levels = c("all", "het", "enough_counts"))]

saveRDS(m_dt, "Output/mae_freqs.Rds")

#' ## Plot
g <- ggplot(m_dt, aes(type, value)) + geom_boxplot() + theme_bw() + 
    labs(y = "SNVs per patient", x = "Filter cascade")
ggplotly(g)

m_dt[value > 1e5]
