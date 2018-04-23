suppressPackageStartupMessages(source("src/r/config.R"))

exons_gene_dt <- readRDS("./resources/exons_gene_dt.Rds")

# Extract first exons
setorder(exons_gene_dt, start)
plus <- exons_gene_dt[strand == "+", .SD[1], by = gene_id]
setorder(exons_gene_dt, -end)
minus <- exons_gene_dt[strand == "-", .SD[1], by = gene_id]
first <- rbind(plus, minus)

columns <- c("seqnames", "start", "end", "gene_id", "exon_id", "strand")
first <- first[seqnames != "MT"]

# write to bed file
write.table(first[, columns, with = F], "resources/sample_first_exon.bed", quote = F, row.names = F, col.names = F)

#se_strand <- readRDS("/s/project/mitoMultiOmics/processed_expression/counts_batch3.Rds")
se <- readRDS("/s/project/genetic_diagnosis/processed_data/Rds/se_batch3.Rds")
#se <- readRDS(file.path(PROC_DATA, "Rds/se_batch3.Rds"))

# get exon counts
se_ex <- assay(se) %>% as.data.table()
se_ex[, exon_id := 1:.N]
se_ex <- merge(se_ex, exons_gene_dt, by = "exon_id")

# remove genes without counts
se_gene <- se_ex[, lapply(.SD, sum), by = gene_name, .SDcols = colnames(assay(se))]
se_gene[, sum_count := rowSums(.SD[, !"gene_name"])]
zero_genes <- se_gene[sum_count == 0, gene_name]
se_ex <- se_ex[!gene_name %in% zero_genes]

sample_ids <- grep("[0-9]", names(se_ex), value = T)
first <- merge(first, se_ex[, c("exon_id", sample_ids), with = F], by = "exon_id")
first[, mean_counts := rowMeans(.SD[, sample_ids, with = F])]

first[, col := "255,0,0"] # dummy color
first

# write to bed file
columns <- c(columns, "start", "end", "col", "mean_counts")
write.table(first[, columns, with = F], file = "resources/sample_first_exon_counts.bed", quote = F, row.names = F, col.names = F)
