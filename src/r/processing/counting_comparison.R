library(ggpval)
library(ggthemes)
ss2 <- read.table(file.path(PROC_DATA, "rna_batch2_strand_specific.txt"), check.names = F)
ss3 <- read.table(file.path(PROC_DATA, "rna_batch3_strand_specific.txt"), check.names = F)
nss_matrix <- read.table(file.path(PROC_DATA, "rna_batch1_2_non_strand_specific.txt"), check.names = F)
# ss_reversed <- readRDS("./resources/rna_count_matrix_strand_specific_reversed.Rds")
ss_new <- read.table(file.path(PROC_DATA, "rna_batch0_new_annotation.txt"), check.names = F)
ss_old <- read.table(file.path(PROC_DATA, "rna_batch0_old_annotation.txt"), sep = "\t", check.names = F)

sa = fread("/data/ouga/home/ag_gagneur/yepez/Desktop/batch_1_2_3.csv")
sa[FIBROBLAST_ID == "", FIBROBLAST_ID := NA]
sa[GENDER == "", GENDER := NA]
sa[EXOME_ID == "", EXOME_ID := NA]

sa0 = SAMPLE_ANNOTATION[RNA_ID %in% colnames(ss0)][, .(EXOME_ID, RNA_ID, FIBROBLAST_ID, GENDER, TISSUE, RNA_PERSON)]
sa0[, BATCH := "B0"]
sa_all = rbind(sa0, sa, use.names = T)
sa_all[TISSUE == "fibroblasts", TISSUE := "FIBROBLAST"]
write.table(sa_all, file.path(PROC_DATA, "sample_annotation_0to3.txt"), quote = F, col.names = T, row.names = F, sep = ",")


batch2 = sa[BATCH == "B2", RNA_ID]
batch1 = sa[BATCH == "B1", RNA_ID]
batch1_fib = sa[BATCH == "B1" & TISSUE == "fibroblasts", RNA_ID]

nss1 = nss_matrix[, batch1_fib]
write.table(nss1, file.path(PROC_DATA, "rna_batch1.txt"), quote = F, col.names = T, row.names = T)
nss2 = nss_matrix[, batch2]

# Cut to include genes whose 95% have more than 10 reads only
filter_percentage <- function(count_matrix, percentage = .95, min_reads = 10){
    count_matrix[apply(count_matrix, 1, quantile, percentage) > min_reads, ]
}

ss95 = filter_percentage(ss2)
dim(ss95) # 16001

nssB195 = filter_percentage(nss1)
dim(ssB195) # 19283 with all tissues, 16297 only fibroblasts

nss95 = filter_percentage(nss2)
dim(nss95)   # 16507

ss_new = filter_percentage(ss_new)
dim(ss_new)  # 15285

ss_old = filter_percentage(ss_old)
dim(ss_old)  # 12606

ss395 = filter_percentage(ss3)
dim(ss395)  # 15436


medians_dt = data.table(gene_name = row.names(ss_new), Median = rowMedians(as.matrix(ss_new)) )
medians_dt = left_join(medians_dt, et, by = "gene_name") %>% as.data.table

png("/data/ouga/home/ag_gagneur/yepez/LRZ Sync+Share/LMU/Group_meetings/2017_12_journal_club/transcript_boxplot.png",
    width = 700, height = 500, res = 140)
ggplot(medians_dt[transcript_type %in% c("protein_coding", "lincRNA", "pseudogene", "antisense")], aes( transcript_type, Median)) + 
    geom_boxplot() + scale_y_log10() + theme_bw(base_size = 14)
dev.off()

et = unique(exons_gene_dt[, .(gene_name, transcript_type)]) 
unique(et$transcript_type)
table(et$transcript_type)

old_genes <- intersect(rownames(ss_old), et[transcript_type != "protein_coding", gene_name]) # 537
new_genes <- intersect(rownames(ss_new), et[transcript_type != "protein_coding", gene_name]) # 3623
table_old <- table(et[gene_name %in% old_genes & transcript_type != "protein_coding", transcript_type]) %>% as.data.table
table_old[, Count_Type := "Old"]
table_new <- table(et[gene_name %in% new_genes & transcript_type != "protein_coding", transcript_type]) %>% as.data.table
table_new[, Count_Type := "New"]
table_both = rbind(table_old, table_new)
setnames(table_both, "V1", "Transcript_type")
table_both[, Count_Type := factor(Count_Type, levels = c("Old", "New"))]

png("/data/ouga/home/ag_gagneur/yepez/LRZ Sync+Share/LMU/Group_meetings/2017_12_journal_club/nc_barplot.png",
    width = 1400, height = 700, res = 140)
ggplot(table_both[N>3], aes(Transcript_type, N, fill = Count_Type)) + 
    geom_bar(stat = "identity", position = "dodge") + 
    geom_text(aes(label = N), color="black",
              position = position_dodge(width = .9), size=4, hjust = -.1) +
    coord_flip() + theme_minimal(base_size = 14) +
    scale_fill_brewer(palette="Paired") + scale_y_log10()
    
dev.off()

hist(log10(nss2+1), breaks=100)
hist(log10(ss95 + 1), breaks=100)

png(file.path(DIR_retreat, "hist_b2_nss.png"), width = 600, height = 500, res = 115)
hist(log10(nss95 + 1), breaks=100, xlab = "B2, Non Strand Specific", cex.axis = 1.4, cex.lab = 1.4)
dev.off()

png(file.path(DIR_retreat, "hist_b2_ss.png"), width = 600, height = 500, res = 115)
hist(log10(ss95 + 1), breaks=100, xlab = "B2, Strand Specific", cex.axis = 1.4, cex.lab = 1.4)
dev.off()

png(file.path(DIR_retreat, "hist_b1_nss.png"), width = 600, height = 500, res = 115)
hist(log10(ssB195 + 1), breaks=100, xlab = "B1, Non Strand Specific", cex.axis = 1.4, cex.lab = 1.4)
dev.off()

# Scatterplot medians
png(file.path(DIR_retreat, "scatterplot_strand_counts.png"), width = 600, height = 600, res = 120)
plot(colMedians(nss95[both_genes, ]), colMedians(ss95[both_genes,]), log = "xy", xlab = "Median Gene Count, Non Strand Specific", ylab = "Median Gene Count, Strand Specific")
abline(0,1)
dev.off()

dt_sums = data.table(mean = c(colMeans2(nss95[both_genes, ]), colMeans2(ss95[both_genes, ])),
                     type = c(rep("Non Strand Specific", ncol(nss95[both_genes,])), rep("Strand Specific", ncol(ss95[both_genes, ])) )
)

png(file.path(DIR_retreat, "boxplot_strand_counts.png"), width = 700, height = 600, res = 120)
g = ggplot(dt_sums, aes(type, mean)) + geom_boxplot() + labs(y = "Mean Gene Count Per Patient") + theme_bw(base_size = 16)
wilcox.test(colMeans2(nss95[both_genes, ]), colMeans2(ss95[both_genes, ]), paired = T)  # 1.421e-14
add_pval(g, pairs = list(c(1,2)), annotation = "P = 1.4e-14")
dev.off()


rownames(ss95) %in% protein_coding_genes %>% sum
rownames(nss95) %in% protein_coding_genes %>% sum


nss_genes = setdiff(rownames(nss95), rownames(ss95))
both_genes = intersect(rownames(nss95), rownames(ss95))
png(file.path(DIR_retreat, "hist_b2_nss_sub.png"), width = 600, height = 500, res = 115)
hist(log10(nss95[nss_genes, ] + 1), breaks=100, xlab = "B2, Non Strand Specific Only", cex.axis = 1.4, cex.lab = 1.4)
dev.off()
