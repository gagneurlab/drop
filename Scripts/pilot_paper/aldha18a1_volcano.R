#'---
#' title: ALDH18A1 with volcano plot
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

#' # FILE 

#' paper data
file_80256_mae <- file.path(TIDYDIR, "ext_suppl_80256_mae_deseq_results.tsv")
tmp <- fread(file_80256_mae)
tmp[hgncid=='ALDH18A1']

#' all stats for all patients
file_mae <- file.path(TIDYDIR, "mae_deseq_wes_stats.RDS")
dt <- readRDS(file_mae)
head(dt)


#' # Plot
#' 

patdt= dt[FIBROBLAST_ID==80256 ]
colnames(patdt)

#' Check how much is filtered
summary(patdt$exacmaf)
summary(patdt$coverage)

#' Assign colors
#' 
patdt[,snv_filter:='expressed']
patdt[ padj<0.05 & (alt_allele_freq>=0.8 | alt_allele_freq<=0.2),
    snv_filter:='MAE'
]
patdt[ snv_filter=='MAE' & (is.na(exacmaf) | exacmaf< 0.001), 
    snv_filter:='rare MAE'
]


#+ message=F
plot_ly(
    data=patdt[snv_filter=='expressed'],
    name='RNA>10',
    marker=list(color='grey'),
    type='scatter', mode='markers',
    x=~alt_allele_freq,
    y=~-log10(pvalue),
    text=~paste(HGNCID, "<br>REF:", ref, ", ALT:", alt, "<br>ExAC_MAF", exacmaf)
)%>% 
    add_trace(
        data=patdt[snv_filter=='MAE'],
        name='MAE',
        marker=list(color='blue')
    ) %>%
    add_trace(
        data=patdt[snv_filter=='rare MAE'],
        name='rare MAE',
        marker=list(color='red')
    ) %>%
    add_annotations(
        data=patdt[HGNCID=='ALDH18A1'],
        ax=~ifelse(alt_allele_freq<0.5, 40, -40), ay=-40
    ) %>%
    layout(
        xaxis=list(title='Alternative allele frequency'),
        yaxis=list(title='-log10( P-value )')
    )





