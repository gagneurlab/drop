#'---
#' title:  
#' author: Daniel Bader
#' wb:
#'   input: 
#'   output: 
#' output: 
#'   html_document:
#'     toc_float: yes
#'     toc: yes
#'---
#'

#+ echo=F
source("src/r/config.R")


#' # Read tidy data
#+ 
file_tidy_rna <- file.path(DATADIR, "processed_expression/rna_aberrant_expression.RDS")
rna_dt <- readRDS(file_tidy_rna) %>% as.data.table()
rna_dt[,FIBROBLAST_ID := as.character(FIBROBLAST_ID)]

#' tidy proteomics
#+
file_tidy_pichler_100min <- file.path(PROC_DATA, "proteome_pichler_100min.tsv")
# prot_dt <- fread(file_tidy_pichler_100min)
# head(prot_dt)

#' ## paper proteomics
#+
file_proteome_bundle <- file.path(DATADIR, "processed_expression/proteome_normalize_fibro_bundle.Rdata")
load(file_proteome_bundle, verbose = T)
head(proteome_limma_res_list[[1]])


#' transform to tidy
#' 
#+
proteome_tidy_dt <- data.table()
for(fib in names(proteome_limma_res_list)){
    dt <- proteome_limma_res_list[[fib]]
    setnames(dt, sub('^gene','hgncid',names(dt)))
    dt[,FIBROBLAST_ID:=as.character(fib)]
    proteome_tidy_dt <- rbind(proteome_tidy_dt, dt)
}
head(proteome_tidy_dt)


#' # Complex 1 genes
#+
extsuppl_prefix= '/s/project/mitoMultiOmics/paper_nature_genetics/paper_supplement_data/ext_suppl_'
FILE_corum_database = paste0(extsuppl_prefix, "corum_mammalian_protein_complexes.tsv")
corumdt = as.data.table(read.delim(FILE_corum_database, stringsAsFactor=F))

genes_complex1 <- corumdt[
            grepl('Respiratory chain complex I', Complex.name), 
            strsplit(gene, ',')[[1]] ] %>% sort %>% unique
#' 
#' # RNA vs Protein fold change plots
#' 

rnaprot_dt <- get_table_for_rna_vs_prot_fc_plot(rna_dt, proteome_tidy_dt)

rnaprot_dt[,plotted_prot_fc:=prot_log2fc]
reliable_freq <- 0.5


#+ 

myfcplot <- function(
    rnaprot_dt,
    fib='73804',
    gene='MGST1',
    fill_na_fc=0.008,
    xpos_label= 0.13
){
    
    #' fill reliable NA
    rnaprot_dt[
        FIBROBLAST_ID==fib & prot_freq_na<reliable_freq & is.na(prot_log2fc), 
        plotted_prot_fc:=log2(fill_na_fc)
        ]
    
    LSD::heatscatter(
        x=2^rnaprot_dt[FIBROBLAST_ID==fib, log2FoldChange],
        y=2^rnaprot_dt[FIBROBLAST_ID==fib, plotted_prot_fc], 
        colpal='blues',
        log='xy',
        xlab='RNA fold change',
        ylab='Protein fold change',
        main=paste0("Patient #",fib),
        las=1,
        cex.lab=1.5, yaxt='n'
    )
    abline(h=0.01); grid()
    axis_at <- c(fill_na_fc, 0.05, 0.1, 0.5, 1, 2, 5, 10)
    axis(side = 2, at = axis_at, labels = sub(fill_na_fc, 'not detected', axis_at), las=1)
    
    tmpdt <- rnaprot_dt[FIBROBLAST_ID==fib & hgncid %in% genes_complex1]
    points(
        x=2^tmpdt[, log2FoldChange],
        y=2^tmpdt[, plotted_prot_fc],
        col='red',
        pch=15
    )
    
    tmpdt <- rnaprot_dt[FIBROBLAST_ID==fib & hgncid==gene]
    text(
        x=2^tmpdt[, log2FoldChange],
        y=2^tmpdt[, plotted_prot_fc], 
        labels = gene, pos=4
    )
}

#+ rnaprot_73804_mgst1, fig.width=8, fig.height=8
par(mar=c(6,6,4,1))
myfcplot(rnaprot_dt, '73804', 'MGST1', 0.008)

#+ rnaprot_80256_aldh18a1, fig.width=8, fig.height=8
par(mar=c(6,6,4,1))
myfcplot(rnaprot_dt, '80256', 'ALDH18A1', 0.008)


#' ## Plotly version
#' 
#' to be extended
#' 
p <- plot_ly(
    data=rnaprot_dt[FIBROBLAST_ID==fib],
    type='scatter',
    mode='markers',
    x=~2^log2FoldChange, 
    y=~2^plotted_prot_fc,
    text=~hgncid
)%>% layout(xaxis=list(type="log"), yaxis=list(type="log"))

#+ warnings=F
p






