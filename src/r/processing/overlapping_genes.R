library(ggpval)
library(ggthemes)
ss_matrix <- read.table(file.path(PROC_DATA, "rna_batch2_strand_specific.txt"))
nss_matrix <- read.table(file.path(PROC_DATA, "rna_batch1_2_non_strand_specific.txt"))


neg = genes_en[strand(genes_en) == "-", ]
strand(neg) <- "*"
pos = genes_en[strand(genes_en) == "+", ]
strand(pos) <- "*"
hits = findOverlaps(neg, pos)
hits

expressed_genes <- row.names(ssr95) %>% unique

unique(gene_mapping$transcript_type)
table(gene_mapping$transcript_type) #20345 protein coding
protein_coding_genes <- gene_mapping[transcript_type == "protein_coding", unique(gene_name)]
non_protein_coding_genes <- gene_mapping[transcript_type != "protein_coding", unique(gene_name)]
miRNAs <- gene_mapping[transcript_type == "miRNA", unique(gene_name)]


setkey(gene_mapping, "gene_id")
gene_names_pos <- gene_mapping[pos[to(hits)]$gene_id, gene_name]
# pos_within_ss_reversed <- gene_names_pos %in% rownames(ss_reversed)

gene_names_neg <- gene_mapping[neg[from(hits)]$gene_id, gene_name]
# neg_within_ss_reversed <- gene_names_neg %in% rownames(ss_reversed)

ss_reversed[gene_names_pos[pos_within_ss_reversed], ]
ss_reversed[gene_names_neg[neg_within_ss_reversed], ]

ssr

res <- data.table(pos_genes = gene_names_pos, neg_genes = gene_names_neg, pos = logical(length(hits)), neg = logical(length(hits)))

nss95 = nss2[apply(nss2, 1, quantile, .95) > 10, ] 
dim(nss95)   # 16507

ssr95 = ssr[apply(ssr, 1, quantile, .95) > 10, ] 
dim(ssr95)   # 16001
subset = intersect(protein_coding_genes, rownames(ssr95))
subset2 = intersect(non_protein_coding_genes, rownames(ssr95))

res[, pos := pos_genes %in% subset, by = .I]
res[, neg := neg_genes %in% subset, by = .I]
res[, pos_nonc := pos_genes %in% subset2, by = .I]
res[, neg_nonc := neg_genes %in% subset2, by = .I]

nrow(res[pos_nonc == T & neg_nonc == T])
resT = res[pos == T & neg == T, ]

dif_genes = setdiff(rownames(nss95), rownames(ssr95))
hist(nss95[dif_genes, ], log = "x")
hist(log10(nss95[dif_genes, ] +1), breaks=100)
hist(log10(nss95 +1), breaks=100)

plot(colMedians(nss95[resT$pos_genes, ]), colMedians(ssr95[resT$pos_genes, ]) ); abline(0,1)
plot(colMedians(nss95[resT$neg_genes, ]), colMedians(ssr95[resT$neg_genes, ]) ); abline(0,1)

# plot(rowSums2(nss95[resT$pos_genes, ]), rowSums2(ssr95[resT$pos_genes, ]), log = "xy" ); abline(0,1)
# plot(rowSums2(nss95[resT$neg_genes, ]), rowSums2(ssr95[resT$neg_genes, ]), log = "xy" ); abline(0,1)




neg_e = exons_en[strand(exons_en) == "-", ]
strand(neg_e) <- "*"
pos_e = exons_en[strand(exons_en) == "+", ]
strand(pos_e) <- "*"
hits_e = findOverlaps(neg_e, pos_e)
hits_e

exon_se_matrix = as.matrix(se_ex[, batch2, with = F])
exon_sn_matrix = as.matrix(se_exn[, batch2, with = F])
exon_ss95 = exon_se_matrix[apply(exon_se_matrix, 1, quantile, .95) > 10, ]
rm(exon_se_matrix)
exon_sn95 = exon_sn_matrix[apply(exon_sn_matrix, 1, quantile, .95) > 10, ] 
rm(exon_sn_matrix)

se2 = se_ex[transcript_type == "protein_coding"]$exon_id
sen2 = se_exn[transcript_type == "protein_coding"]


exon_names_pos <- pos_e[to(hits_e)]$exon_id
exon_names_neg <- neg_e[from(hits_e)]$exon_id
res_e <- data.table(pos_exons = exon_names_pos, neg_exons = exon_names_neg, pos = logical(length(hits_e)), neg = logical(length(hits_e)))



