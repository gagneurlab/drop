#'---
#' title: Kuester data QC
#' author: Daniel Bader
#' wb:
#'   input: [ 
#'     "/s/project/mitoMultiOmics/raw_data/proteome/20170614_kopajtich_kuester_proteome/m1_lfq_single_shot_proteinGroups.txt",
#'     "/s/project/mitoMultiOmics/raw_data/proteome/20170614_kopajtich_kuester_proteome/m2_lfq_rphp_proteinGroups.txt",
#'     "/s/project/mitoMultiOmics/raw_data/proteome/20170614_kopajtich_kuester_proteome/m3_lfq_tmt_proteinGroups.txt",
#'     "/s/project/mitoMultiOmics/raw_data/proteome/20170614_kopajtich_kuester_proteome/m4_lfq_id_trinity_proteinGroups.txt"
#'     ]
#'   output: "/s/project/genetic_diagnosis/processed_data/proteome_kuester_method_trial_lfq.tsv"
#' output: 
#'   html_document:
#'     toc_float: yes
#'     toc: yes
#'---
#'

#+ echo=F
source("src/r/config.R")
#+ 
opts_chunk$set(message=T)

#' # Data
#' 
#' Kuester proteome files:
protdir <- file.path(RAWDIR, "proteome", "20170614_kopajtich_kuester_proteome")
files_kuester <- list.files(protdir, pattern = "^m.*txt$", full.names = T)
print(files_kuester)

file_out <- file.path(PROC_DATA, "proteome_kuester_method_trial_lfq.tsv")


#' 
#' # Single Shot
#' 
tmp_file_kuester <- files_kuester[1]

#' 
#' ## LFQ
#'
#' * extract intensities from proteinGroups.txt
#' * create tidy proteome data table "pdt"
#' * filter protein expression by gene properties
#' * compute NA frequencies by gene, sample
#' 

pdt <- wrapper_proteinGroupsTxt_to_tidy_table(tmp_file_kuester)
print(pdt)
pdt[PROTEOME_ID=='nhdf.p9.rep1', PROTEOME_ID:="nhdf.p9"]
pdt[,ms_method:='single_shot']


#' # Read remaining Kuester files
#' 

#' ## rphp
pdt2 <- wrapper_proteinGroupsTxt_to_tidy_table(
    files_kuester[2],
    intensity_column_pattern = "LFQ.intensity."
)
pdt2[PROTEOME_ID=='nhdf.p9.rep1', PROTEOME_ID:="nhdf.p9"]
pdt2[,ms_method:='rphp']
print(pdt2)


#' ## TMT   
pdt3 <- wrapper_proteinGroupsTxt_to_tidy_table(
    files_kuester[3],
    intensity_column_pattern = "Reporter.intensity.corrected."
)
pdt3[,ms_method:='tmt']

#' remove "TMT_10plex_Test_Trinity" samples
pdt3 <- pdt3[!grepl("TMT_10plex_Test_Trinity",PROTEOME_ID)]

#' add sample IDs from separate table
file_tmt_sample_map <- file.path(
    'resources', '201706_kuester_tmt_sample_mapping.tsv'
)
tmt_map <- fread(file_tmt_sample_map)
for(i in 0:9){
    pdt3[PROTEOME_ID==i, PROTEOME_ID:=tmt_map[ID==i, PROTEOME_ID]]
}

print(pdt3)


#' ## Trinity
pdt4 <- wrapper_proteinGroupsTxt_to_tidy_table(
    files_kuester[4],
    intensity_column_pattern = "Intensity."
)
pdt4[, PROTEOME_ID:="nhdf.p9"]
pdt4[,ms_method:='trinity']
print(pdt4)


#' 
#' # Combine all 4 data sets
#' 


#+ echo=F
opts_chunk$set(message=F, echo=F)

#+
pdtall <- rbindlist(list(pdt, pdt2, pdt3, pdt4)) 
write_tsv(pdtall, file = file_out)



#' 
#' ## Sample stats
#' 

ms_summary <- unique(pdtall[,
        .(.N, 
            NA_FREQ_BY_PROTEOME_ID, 
            min_na_freq_by_gene=min(NA_FREQ_BY_GENE_NAME, na.rm=T),
            median_na_freq_by_gene=median(NA_FREQ_BY_GENE_NAME, na.rm=T),
            mean_na_freq_by_gene=mean(NA_FREQ_BY_GENE_NAME, na.rm=T),
            max_na_freq_by_gene=max(NA_FREQ_BY_GENE_NAME, na.rm=T),
            min_lfq=min(LFQ_INTENSITY, na.rm=T), 
            median_lfq=median(LFQ_INTENSITY, na.rm=T), 
            mean_lfq=mean(LFQ_INTENSITY, na.rm=T), 
            max_lfq=max(LFQ_INTENSITY, na.rm=T)
        ), 
        by=c("ms_method", "PROTEOME_ID")
    ]
)


#+
DT::datatable(ms_summary, filter='top', rownames = F)

#+
plot_ly(data=unique(ms_summary[,.(ms_method, N)]),
    type='bar',
    x=~ms_method,
    y=~N,
    color=~ms_method,
    showlegend=F
)

#+
plot_ly(data=ms_summary,
        type='bar',
        x=~paste(ms_method,PROTEOME_ID),
        y=~NA_FREQ_BY_PROTEOME_ID,
        color=~ms_method
)%>% layout(
    xaxis=list(categoryorder='array', categoryarray=~ms_method, title=''),
    margin=list(b=100)
)

#'
#' ## Intensity distribution
#' 


#+
plot_ly(data=pdtall,
        type='box',
        x=~paste(ms_method,PROTEOME_ID),
        y=~LFQ_INTENSITY,
        color=~ms_method
)%>% layout(
    xaxis=list(categoryorder='array', categoryarray=~ms_method, title=''),
    yaxis=list(type='log',exponentformat='power'),
    margin=list(b=100)
)

#+
p <- ggplot(
        pdtall, 
        aes(x=paste(ms_method,PROTEOME_ID), y=LFQ_INTENSITY)
    ) + 
    geom_violin() + 
    scale_y_log10() + 
    labs(x='') + theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplotly(p)


#' ## Investigate Dynamic range
#+
plot_ly(
        data= pdtall[,
            .(inter_quartile_range_log10 = IQR(x=log10(LFQ_INTENSITY), na.rm=T)), 
            by=ms_method
        ],
        x=~ms_method,
        y=~inter_quartile_range_log10
    )%>% 
    add_bars(
        color=~ms_method
    )

#+
plot_ly(
        data= pdtall[, .(q90_range_log10 = I(
                    quantile(x=log10(LFQ_INTENSITY), probs=0.95, na.rm=T) -
                    quantile(x=log10(LFQ_INTENSITY), probs=0.05, na.rm=T)
                )
            ), 
            by=ms_method
        ],
        x=~ms_method,
        y=~q90_range_log10
    )%>% 
    add_bars(
        color=~ms_method
    )



#'
#' ## Missing values per detected protein


#+
plot_ly(
        data = pdtall[,
            .(mean_na_freq_by_gene=mean(NA_FREQ_BY_GENE_NAME)), by=ms_method
        ],
        x=~ms_method,
        y=~mean_na_freq_by_gene
    )%>% 
    add_bars(
        color=~ms_method
    )


#+
p <- plot_ly(
        data=pdtall, 
        x=~paste(ms_method,PROTEOME_ID),
        y=~NA_FREQ_BY_GENE_NAME
    ) %>%
    add_boxplot(
        color=~ms_method, boxpoints=FALSE
    )%>% layout(
        xaxis=list(categoryorder='array', categoryarray=~ms_method, title=''),
        margin=list(b=100)
)
p

#+
p <- ggplot(
        pdtall, 
        aes(x=paste(ms_method,PROTEOME_ID), y=NA_FREQ_BY_GENE_NAME)
    ) + 
    geom_violin() + 
    labs(x='') + theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplotly(p)



#' 
#' # Compare NHDF betw methods
#' 

#' ## Heatmap
#' 
#+
pdt_nhdf <- dcast(
    pdtall[PROTEOME_ID=='nhdf.p9'], 
    GENE_NAME ~ ms_method, 
    value.var = 'LFQ_INTENSITY'
)
pmat_nhdf <- as.data.frame(pdt_nhdf)
rownames(pmat_nhdf) <- pmat_nhdf[,'GENE_NAME']
pmat_nhdf <- pmat_nhdf[,setdiff(colnames(pmat_nhdf),'GENE_NAME')]

#+
tmp= cor(pmat_nhdf, use="pairwise.complete.obs")
plot_ly(z=tmp, type='heatmap', 
    x=rownames(tmp), y=colnames(tmp),
    colors=c('grey95', 'dodgerblue')
)%>% layout(title="Pearson correlation")


#+
tmp= cor(pmat_nhdf, use="pairwise.complete.obs", method='spear')
plot_ly(
        z=tmp, type='heatmap', 
        x=rownames(tmp), y=colnames(tmp),
        colors=c('grey95', 'dodgerblue')
    )%>% layout(
        title="Spearman rank correlation"
    )

#' 
#' ## Scatter protein vs protein
#' 
#+
mydiag <- list(
    type = "line",
    line = list(color = "pink"),
    xref = "x",
    yref = "y",
    x0 = 1e4, y0 = 1e4,
    x1 = 1e12, y1 = 1e12
)
nhdf_trinity <- plot_ly(
        data=pdt_nhdf, x=~trinity, text=~GENE_NAME, alpha=0.2
    ) %>%
    layout(
        shapes=list(mydiag),
        xaxis=list(type='log'),
        yaxis=list(type='log')
    )

#+ rphp vs trinity
p1 <- nhdf_trinity %>% add_markers(y=~rphp)
p1

#+
plot(rphp ~ trinity, data=pdt_nhdf, log='xy')
abline(0,1)


#+ tmt vs trinity
nhdf_trinity <- plot_ly(
    data=pdt_nhdf, x=~trinity, text=~GENE_NAME, alpha=0.2
) %>%
    layout(
        shapes=list(mydiag),
        xaxis=list(type='log'),
        yaxis=list(type='log')
    )
nhdf_trinity %>% add_markers(y=~tmt)


#+
x=as.matrix(na.omit(log2(pdt_nhdf[,c(5,4), with=F])))
pca=prcomp(x)
a= mean(x[,2], na.rm=T)- (0.51/0.86)*mean(x[,1], na.rm=T)

plot(x, xlim=c(15,40), ylim=c(5,30))
abline(0,1)
abline(a=a , b=0.51/0.86)



#' 
#' # Correlation within TMT data
#' 

#' * Spearman rank correlation with use=="complete.obs" 
#'   then missing values are handled by casewise deletion 
#'   (and if there are no complete cases, that gives an error). 
#' * Hierarchical clustering.

tmt_log2fc_intensity = normalize_expression_matrix(
    convert_tidy_table_to_numeric_matrix(
        tidytable = pdtall[ms_method=='tmt'], 
        coln_rownames = 'GENE_NAME',
        coln_colnames = 'PROTEOME_ID',
        coln_values = 'LFQ_INTENSITY'
    ), 
    sizefactor = T, rowcenter = T, log2scale = T, nonzero = 0
)
tmt_cormat <- cor(tmt_log2fc_intensity, method = 's', use = 'complete.obs')

#+
# hc <- as.dendrogram(hclust(dist(tmt_cormat)))
# heatmap_notrace(tmt_cormat, Rowv=hc, Colv=hc)

#+ fig.width=8, fig.height=8
heatmaply_cor(tmt_cormat)

#+
# END
