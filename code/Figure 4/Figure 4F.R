mat=read.table("DML.CpG.TSS.deta20.NAfil.matrix.txt")
group_list<-rep(c("Adjacent","Tumor"),13)
names(group_list)=colnames(mat)
annotation_col=data.frame(group_list)

mtdata<-matrix(as.numeric(unlist(mat)),nrow = length(rownames(mat)),ncol =length(colnames(mat)))
rownames(mtdata)=rownames(mat)
colnames(mtdata)=colnames(mat)
a=rowSums(mtdata)
mtdata2=mtdata[-which(a==0),]
library("RColorBrewer")
colset <- brewer.pal(5, "Set1")


anno_col=list(group_list=c(Tumor=colset[1],Adjacent=colset[2]))

pheatmap(mtdata2, cluster_cols = TRUE,
annotation_col=annotation_col,clustering_method = "ward.D",
show_colnames = F,show_rownames=F,scale="row",
annotation_colors=anno_col,
color = colorRampPalette(c("#3366CC", "white", "darkred"))(250))
