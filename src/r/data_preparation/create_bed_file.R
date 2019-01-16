
library(GenomicAlignments)

myRange <- GRanges(
    "chr17",
    IRanges(73894519, 73901314)
)
myRange

myCounts <- summarizeOverlaps(myRange, BamFile("/s/project/mitoMultiOmics/raw_data/helmholtz/MUC3202/RNAout/paired-endout/stdFilenames/MUC3202.bam.bai"))
myCounts


red <- reduce(gencode_op)
red
red_dv <- unlist(red)
red_dt <- as.data.table(red_dv)
red_dt$gene_id = names(red_dv)
red_dt <- red_dt[order(seqnames, start, end)]
cat(red_dt[,paste0(seqnames, "\t", start, "\t", end, "\t", gene_id)], sep="\n", file = "resources/reduced_v29.bed")
unlink("resources/reduced_v29.bed.idx")

fread("resources/reduced_v29.bed")

