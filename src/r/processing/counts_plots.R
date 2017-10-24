library("pheatmap")
library("RColorBrewer")

DIR_retreat = "/data/ouga/home/ag_gagneur/yepez/LRZ Sync+Share/LMU/Group_meetings/2017_09_lab_retreat"

sa = SAMPLE_ANNOTATION[RNA_ID %in% samples]
sa = fread("/data/ouga/home/ag_gagneur/yepez/Desktop/batch_1_2.csv")
sa[, BATCH := as.factor(BATCH)]
sa_batch = data.frame(batch = sa$BATCH)
row.names(sa_batch) = sa$RNA_ID

sa[, TISSUE := as.factor(TISSUE)]
sa_tissue = data.frame(tissue = sa$TISSUE)
row.names(sa_tissue) = sa$RNA_ID


ann_colors <- list(batch = c(B1 = "brown", B2 = "darkgreen"),
                   tissue = c(fibroblasts = "sienna1", heart = "firebrick1", kidney = "blueviolet", liver = "darkgray", muscle = "gold") )

# Heatmap of distances
sampleDists <- dist(t(log10(se_matrix_n + 1)))
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

pheatmap(as.matrix(sampleDists),
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors,
         border_color = NA,
         annotation_row = sa_batch,
         annotation_colors = ann_colors,
         annotation_col = sa_tissue,
         # annotation_colors = tissue_colors,
         show_colnames = F,
         fontsize = 16, fontsize_row = 10,
         main = "All samples", 
          filename = file.path(DIR_retreat, "heatmap_batches.png")
         )


