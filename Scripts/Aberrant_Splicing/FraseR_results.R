#'---
#' title: FraseR results up to batch 5
#' author: Christian Mertes
#' wb:
#'  input:
#'   - fds: "/s/project/fraser/analysis/datasets/savedObjects/prokisch_batch5/fds-object.RDS"
#'---
suppressMessages({
    library(ggplot2)
    devtools::load_all("../FraseR")
    # devtools::load_all("../fraser")
})
source("../genetic_diagnosis/Scripts/_functions/gene_annotation/add_gene_info_cols.R")

fdsFile <- "/s/project/fraser/analysis/datasets/savedObjects/prokisch_batch5/fds-object.RDS"
fdsFile <- snakemake@input$fds

fdsPath <- dirname(dirname(dirname(fdsFile)))
fdsName <- basename(dirname(fdsFile))


#'
#' # FraseR results
#'
#' This contains the prokisch dataset up to batch 5
#'
fds <- loadFraseRDataSet(fdsPath, fdsName)

#'
#' ## Sample annotation
#'
DT::datatable(as.data.table(colData(fds)))

#'
#' ## Results by sample and gene
#'
results <- results(fds, fdrCut=0.1, dPsiCut=0.1)


# rename columns and cleaning data, echo=FALSE, results="hide"
colnames(mcols(results))[12] <- "totalExpression"
results$p.adj <- abs(results$p.adj)
results$pvalue <- abs(results$pvalue)

resultsdt <- as.data.table(results)
setnames(resultsdt, "hgnc_symbol", "gene_name")
resultsdt <- add_all_gene_info(resultsdt)

sa <- fread(snakemake@config$SAMPLE_ANNOTATION)
resultsdt[, gene_name := toupper(gene_name)]
resultsdt <- left_join(resultsdt, sa[, .(RNA_ID, FIBROBLAST_ID, EXOME_ID, PEDIGREE, KNOWN_MUTATION,
                             CANDIDATE_GENE, BATCH, COMMENT)],
                 by = c("sampleID" = "RNA_ID")) %>% as.data.table



#'
#' ### Number of Events (genes) by sample
#'
filtdt <- resultsdt[abs(p.adj) <= 0.1 &
        abs(deltaPsi) > 0.1 &
        (totalExpression >= 5) &
        meanTotalCts < 100 * totalExpression]


#'
#' Without intron retention (grouped by gene names)
#'
numEventsBySample <- filtdt[type != "psiSite", .N,by=c("sampleID", "gene_name")][
        , .N, by="sampleID"][, .(rank=rank(N), N, sampleID)][order(rank)]

library(ggpubr)
ggplot(numEventsBySample, aes(x=1:length(N), y=N)) +
    geom_histogram(stat="identity") +
    scale_y_log10() +
    ylab("Number of \naberrantly spliced genes") +
    xlab("Ranked subjects") +
    grids() +
    geom_hline(aes(yintercept=median(N)), col="firebrick") +
    annotate("text", label=paste("Median: ", round(median(numEventsBySample$N))), x=10, y=median(numEventsBySample$N)*1.1)


#'
#' Including intron retention (grouped by gene names)
#'
numEventsBySample <- filtdt[, .N,by=c("sampleID", "gene_name")][
    , .N, by="sampleID"][, .(rank=rank(N), N, sampleID)][order(rank)]
ggplot(numEventsBySample, aes(x=1:length(N), y=N)) +
    geom_histogram(stat="identity") +
    scale_y_log10() +
    ylab("Number of \naberrantly spliced genes") +
    xlab("Ranked subjects") +
    grids() +
    geom_hline(aes(yintercept=median(N)), col="firebrick") +
    annotate("text", label=paste("Median: ", round(median(numEventsBySample$N))), x=10, y=median(numEventsBySample$N)*1.1)


#'
#' ### Histogram of the deltaPSI distribution
#'
hist(filtdt[type != "psiSite", deltaPsi], breaks=100, main="Histogram of significant\ndeltaPSI values")


#'
#' ### Global QQ plot of the P-values for PSI_5
#'
myPvals <- as.matrix(pVals(fds, type="psi5"))
myPvalsdt <- data.table(e=ppoints(prod(dim(myPvals))), o=sort(as.vector(abs(myPvals))))
ggplot(myPvalsdt[o < 0.7], aes(x=-log10(e), y=-log10(o))) +
    geom_hex() +
    geom_abline(slope=1, intercept=0, col="firebrick") +
    grids()


# define result functions, echo=FALSE, results="hide"
group_res_by_genes <- function(dt, type=c("psi5", "psi3", "psiSite")){
    local_type <- type
    dt <- dt[type %in% local_type]
    dt <- dt[order(p.adj)]

    dt[,c("rank_by_gene_sample", "num_events_by_gene"):=list(1:.N, .N), by="sampleID,gene_name"]
    dt <-dt[rank_by_gene_sample == 1]
    dt[,rank_by_gene_sample:=NULL]
    return(dt)
}


#'
#' Final results without intron retention (grouped by gene names)
#'
DT::datatable(group_res_by_genes(filtdt, c("psi3", "psi5")))



#'
#' Final results with intron retention (grouped by gene names)
#'
DT::datatable(group_res_by_genes(filtdt, c("psi3", "psi5", "psiSite")))


#' # Full results table for download
write.table(filtdt, "/s/public_webshare/project/genetic_diagnosis/results/FraseR_full_results.tsv", row.names=FALSE, sep="\t")
write.table(group_res_by_genes(filtdt, c("psi3", "psi5")), "/s/public_webshare/project/genetic_diagnosis/results/FraseR_psi3_5_by_gene_results.tsv", row.names=FALSE, sep="\t")
write.table(group_res_by_genes(filtdt, c("psi3", "psi5", "psiSite")), "/s/public_webshare/project/genetic_diagnosis/results/FraseR_psi3_5_site_by_gene_results.tsv", row.names=FALSE, sep="\t")

#' [Download FraseR full results table](https://i12g-gagneurweb.informatik.tu-muenchen.de/project/genetic_diagnosis/results/FraseR_full_results.tsv)
#' 
#' [Download FraseR grouped by genes results table (psi3 and psi5 only)](https://i12g-gagneurweb.informatik.tu-muenchen.de/project/genetic_diagnosis/results/FraseR_psi3_5_by_gene_results.tsv)
#' 
#' [Download FraseR grouped by genes results table](https://i12g-gagneurweb.informatik.tu-muenchen.de/project/genetic_diagnosis/results/FraseR_psi3_5_site_by_gene_results.tsv)


assayNames(fds)
mycor <- cor(as.matrix(assay(fds, "psi5")), method="pear")
pheatmap::pheatmap(mycor)
