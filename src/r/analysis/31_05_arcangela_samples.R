dis_genes <- strsplit(SAMPLE_ANNOTATION[DISEASE == "MITO", unique(KNOWN_MUTATION)], ";") %>% 
    unlist %>% unique %>% na.omit %>% sort


setdiff(disgene_dt[DISEASE == "MITO", HGNC_GENE_NAME], dis_genes)

dis_genes <- union(dis_genes, disgene_dt[DISEASE == "MITO", HGNC_GENE_NAME])


paper_dt <- fread(file.path(dir_gene_info, 'kremer_bader_2016_biorxiv_mitochondrial_disease_genes.tsv'), na.strings=c("NA", ''))
setnames(paper_dt, c('HGNC_GENE_NAME', 'FULL_GENE_NAME', 'MIM_NUMBER', 'OMIM_LINK', 'ENTREZ_GENE_ID'))
paper_dt[, HGNC_GENE_NAME := toupper(HGNC_GENE_NAME)]

setdiff(paper_dt[, HGNC_GENE_NAME], dis_genes)


dis_genes <- union(dis_genes, paper_dt[, HGNC_GENE_NAME])

write.table(dis_genes, "/s/project/mitoMultiOmics/ncRNA/mito_disease_genes.tsv", quote = F, row.names = F, sep = "\n", col.names = F)
length(dis_genes)


library(OUTRIDER)
library(data.table)
library(magrittr)
library(ggplot2)

out_res <- fread("/s/project/scared/results/res_outr_no_filt_prok_allbatches_p_0.05.tsv")
out_res[, V1 := NULL]
out_res_all <- fread("/s/project/scared/results/res_outr_no_filt_allbatches_all.tsv")
out_res_all[, V1 := NULL]
out_res_all[, signif := abs(zScore) > 3 & padjust < .05]
samples <- c("MUC1393", "MUC1395", "MUC2223")
fibs <- c("73638", "73641", "95595")
out_res[sampleID %in% samples]

mat <- out_res_all[, .(sampleID, geneID, rawcounts)]
count_dt <- dcast.data.table(mat, sampleID ~ geneID)
count_mat <- as.matrix(count_dt[,2:ncol(count_dt)])
rownames(count_mat) <- count_dt$sampleID
count_mat[1:5, 1:5]

ods <- OutriderDataSet(countData = count_mat)


out_res_all[sampleID %in% samples][padjust<.99999]
uniqueN(out_res_all$sampleID)

rna_daniel <- readRDS("/s/project/mitoMultiOmics/processed_expression/rna_aberrant_expression.RDS")

fibs %in% rna_daniel$FIBROBLAST_ID
sort(unique(rna_daniel$FIBROBLAST_ID))


rna_daniel[FIBROBLAST_ID %in% fibs][rna_is_signi == T]
?plotVolcano

library(ggthemes)
Ngenes <- uniqueN(out_res_all$geneID)
myvolcano <- function(DT, sample){
    ggplot(DT[sampleID == sample], aes(zScore, -log10(pValue))) + geom_point(aes(col = signif)) +
        theme_bw(base_size = 14) + geom_hline(yintercept = -log10(.05/Ngenes), linetype = "dashed") +
        geom_vline(xintercept = c(-3,3), linetype = "dashed") + ggtitle(sample) + scale_color_ptol()
    
}

cairo_pdf("../../Documents/volcano_plots.pdf", width = 7, height = 5, onefile = T)
myvolcano(out_res_all, "MUC2223")
myvolcano(out_res_all, "MUC1393")
myvolcano(out_res_all, "MUC1395")
dev.off()


g <- ggplot(out_res_all[geneID == "ENSG00000006530.11"], aes(reorder(sampleID, normcounts), normcounts)) + 
    geom_point(aes(col = signif)) + theme_bw(base_size = 14) + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                                                                axis.line = element_line(colour = "black"),
                                                                     axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    labs(x = "Rank", y = "Normalized Counts") + scale_color_ptol() + ggtitle("ENSG00000006530.11")# + scale_y_log10()

png("../../Documents/gene_AGK.png", width = 700, height = 500, res = 120)
g
dev.off()

samples <- c("MUC1393", "MUC1395", "MUC2223")
out_res[sampleID %in% samples]
out_res[geneID == "ENSG00000006530.11"]


g + geom_bar() 



