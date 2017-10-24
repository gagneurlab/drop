library(ggpval)
library(ggthemes)
ss2 <- read.table(file.path(PROC_DATA, "rna_batch2_strand_specific.txt"))
nss_matrix <- read.table(file.path(PROC_DATA, "rna_batch1_2_non_strand_specific.txt"))
# ss_reversed <- readRDS("./resources/rna_count_matrix_strand_specific_reversed.Rds")
sso <- read.table(file.path(PROC_DATA, "rna_batch0_new_annotation.txt"))
ssf <- read.table(file.path(PROC_DATA, "rna_batch0_old_annotation.txt"), sep = "\t", check.names = F)

sa = fread("/data/ouga/home/ag_gagneur/yepez/Desktop/batch_1_2.csv")
batch2 = sa[BATCH == "B2", RNA_ID]
batch1 = sa[BATCH == "B1", RNA_ID]
batch1_fib = sa[BATCH == "B1" & TISSUE == "fibroblasts", RNA_ID]

nss1 = nss_matrix[, batch1_fib]
nss2 = nss_matrix[, batch2]

# Cut to include genes whose 95% have more than 10 reads only
ss95 = ss2[ apply(ss2, 1, quantile, .95) > 10, ] 
dim(ss95) # 16001

nssB195 = nss1[ apply(nss1, 1, quantile, .95) > 10, ] 
dim(ssB195) # 19283 with all tissues, 16297 only fibroblasts

nss95 = nss2[apply(nss2, 1, quantile, .95) > 10, ] 
dim(nss95)   # 16507

sso95 = sso[apply(sso, 1, quantile, .95) > 10, ] 
dim(sso95)  # 16645


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
