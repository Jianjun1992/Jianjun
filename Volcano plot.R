setwd("F:/")
library(Seurat)
library(dplyr)
library(ggplot2)
volcano_data<- read.table("volcano.txt",sep="\t",check.names=F,header=T)
threshold=1
loc_up <-intersect(which(volcano_data$p_val_adj<0.05),which(volcano_data$avg_log2FC>threshold))
loc_down <- intersect(which(volcano_data$p_val_adj<0.05),which(volcano_data$avg_log2FC<(-1*threshold)))
significant <-rep("normal", times=nrow(volcano_data))
significant[loc_up] <- "up"
significant[loc_down] <- "down"
significant <- factor(significant, levels = c("up", "down", "normal"))


plot(x=volcano_data$avg_log2FC, y=-log10(volcano_data$p_val_adj), xlab = "log2(FC)", ylab = "-log10(FDR)", pch=20,size = I(0.2),colour=significant)
i=1
for (i in 1:length(significant)){
  
  if (significant[i]!="normal")
    text(volcano_data$avg_log2FC[i], -log10(volcano_data$p_val_adj[i]), rownames(volcano_data[i,0]),cex = .75,pos=2)
  i=i+1
}


p <-qplot(x=volcano_data$avg_log2FC, y=-log10(volcano_data$p_val_adj), xlab = "log2(FC)", ylab = "-log10(FDR)", size = I(1.5),colour=significant)
p <- p+ scale_color_manual(values = c("up" = "red", "normal" = "black", "down" = "green"))
xline=c(-log2(2),log2(2))
p <- p+geom_vline(xintercept = xline, lty=0,size=I(0.2),colour="grey11")
yline=-log(0.05,10)
p<- p+geom_hline(yintercept = yline, lty=2,size=I(0.2),colour = "grey11")
p<- p+theme_bw()+theme(panel.background = element_rect(colour = "black", size = 1, fill = "white"),panel.grid = element_blank())
png("1.png", width = 3033, height = 2500, res = 500)
print(p)
dev.off()
