library(DESeq2)
setwd("F:/")
mycounts <- read.table("PCA.txt", header = T, row.names = 1)
head(mycounts)
condition <- factor(c(rep("PBS 1d",3), rep("LPS 1d", 3), rep("LPS 7d", 3)), levels = c("PBS 1d","LPS 1d","LPS 7d"))
condition
colData <- data.frame(row.names = colnames(mycounts), condition)
colData
dds <- DESeqDataSetFromMatrix(mycounts, colData, design = ~condition)
dds <- estimateSizeFactors(dds)
raw <- SummarizedExperiment(counts(dds, normalized=FALSE),
                            colData=colData(dds))
nor <- SummarizedExperiment(counts(dds, normalized=TRUE),
                            colData=colData(dds))
vsd <- vst(dds)
rld <- rlog(dds)

plotPCA(vsd)
plotPCA(rld)
dev.off()