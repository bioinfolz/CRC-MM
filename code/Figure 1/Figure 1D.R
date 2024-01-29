data=read.delim("meta.path.res.anno.sort.txt",head=T)
top=data[which(data$I2<=50),]
top=top[which(top$fdr<=0.05),]
top=top[seq(1,20),]

cnt=read.delim("sample.count.txt",head=T,row.names=1,check.names=FALSE)
cnt=cnt[,-1]
sub=cnt[which(rownames(cnt) %in% top$Row.names),]
samp =read.delim("meta.txt",head=T,row.names=2)
case=sub[,which(colnames(sub) %in% rownames(samp[which(samp$study_condition=="CRC"),]))]
control=sub[,which(colnames(sub) %in% rownames(samp[which(samp$study_condition=="control"),]))]
df=cbind(control,case)
des=top[which(top$Row.names %in% rownames(df)),c(1,2)]
rownames(des)=des$Row.names
des=des[rownames(df),]
rownames(df)=des$des
anno=samp[,c(1,2)]
anno=anno[colnames(df),]

library("RColorBrewer")
colset <- brewer.pal(3, "Set1")
colset2 <- brewer.pal(10, "Paired")
ann_colors=list(Study_condition=c(control=colset[1],case=colset[2]),dataset_name=colset2)
pheatmap(df,scale = "row", cluster_cols = FALSE,show_colnames = F,cutree_col=2,
color = colorRampPalette(c("blue", "white", "red"))(255),gaps_col=c(662),
annotation_col=anno,clustering_method = "ward.D")
