setwd("F:/")
library(pheatmap)
library(vegan)
data<-read.table("pheatmap.txt",header=T,sep="\t",row.names=1)
data <- as.data.frame(data) 
data.1 <- decostand(data,"standardize",MARGIN = 1)
View(data.1)
apply(data.1,1,sd)
pheatmap(data.1)