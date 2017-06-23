#'---
#' title: Kuester data QC
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
#+ 
opts_chunk$set(message=T)

#' # Data

#protdir <- file.path(RAWDIR, "proteome", "20170614_kopajtich_kuester_proteome")
protdir <- file.path("~/Documents/tmp_kuester_proteome")

files_kuester <- list.files(protdir, pattern = "^m.*txt$", full.names = T)

#' Kuester proteome files:
print(files_kuester)



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
pdt[,ms_method:='single_shot']


#' 
#' ## iBAQ
#'

#' Same as for LFQ

protein_table_ibaq <- wrapper_proteinGroupsTxt_to_tidy_table(
    tmp_file_kuester, 
    intensity_column_pattern = "iBAQ.", 
    column_intensity = "iBAQ_INTENSITY"
)

print(head(protein_table_ibaq))



#' # Read remaining Kuester files
#' 

#' ## rphp
pdt2 <- wrapper_proteinGroupsTxt_to_tidy_table(
    files_kuester[2],
    intensity_column_pattern = "Intensity.TMT_10plex_Test_"
)
pdt2[,ms_method:='rphp']
print(pdt2)


#' ## TMT   
pdt3 <- wrapper_proteinGroupsTxt_to_tidy_table(
    files_kuester[3],
    intensity_column_pattern = "Intensity."
)
pdt3[,ms_method:='tmt']
print(pdt3)


#' ## Trinity
pdt4 <- wrapper_proteinGroupsTxt_to_tidy_table(
    files_kuester[4],
    intensity_column_pattern = "Intensity."
)
pdt4[,ms_method:='trinity']
print(pdt4)


#' 
#' # Combine all 4 data sets
#' 

#+ echo=F
opts_chunk$set(message=F, echo=F)

#+
pdtall <- rbindlist(list(pdt, pdt2, pdt3, pdt4)) 

#' ## Sample stats

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
plot_ly(data=ms_summary,
    type='bar',
    x=~paste(ms_method,PROTEOME_ID),
    y=~N,
    color=~ms_method
)%>% layout(
    xaxis=list(categoryorder='array', categoryarray=~ms_method, title=''),
    margin=list(b=100)
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
#' ## Gene stats
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

#+
p <- ggplot(
        pdtall, 
        aes(x=paste(ms_method,PROTEOME_ID), y=NA_FREQ_BY_GENE_NAME)
    ) + 
    geom_violin() + 
    labs(x='') + theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplotly(p)

#+
# END
