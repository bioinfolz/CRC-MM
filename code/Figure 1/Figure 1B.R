###pheatmap
library(pheatmap)
setwd("E://13-CRC/submission/revise1/figure 1")
data=read.delim("species.txt",row.names=1)

res=read.delim("merged_metaphlan_table.4lefse.species.txt.res",head=F)
sig=res[which(res$V4>2),]
mat<-data[which(rownames(data) %in% sig$V1),]
log=log(mat+0.000001)
Type = rep(c("Control","CRC"),each=18)
names(Type)=colnames(mat)
annotation_col=data.frame(Type)
library("RColorBrewer")
colset <- brewer.pal(3, "Set1")
ann_colors = list(
Type = c( "Control"= colset[2], "CRC"= colset[1]))

pheatmap(log,scale = "row", cluster_cols = FALSE, 
annotation_col=annotation_col,cutree_rows=2,gaps_col = c(18),
show_colnames = FALSE,annotation_colors=ann_colors)




