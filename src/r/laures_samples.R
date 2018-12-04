library(OUTRIDER)

data <- read.delim2("~/../mertes/TEST_OUTRIDER_COUNTS.tsv")

ods <- OutriderDataSet(countData=data)

ods <- ods[rowMeans(counts(ods)) > 10]
ods
ods <- OUTRIDER(ods, impl="pca")
plotCountCorHeatmap(ods, norm=FALSE)
plotCountCorHeatmap(ods, norm=TRUE)

plotExpressionRank(ods, which(grepl("ENSG00000146282", rownames(ods))), norm=TRUE)
plotExpressionRank(ods, which(grepl("ENSG00000146282", rownames(ods))), norm=FALSE)

sessionInfo()

