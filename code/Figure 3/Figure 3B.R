meta=read.delim("humann_pathabundance.cpm.sum.fil.tsv",row.names=1)
aplicon=read.delim("16s_path_abun_unstrat_descrip.tsv",row.names=1)
mn=sapply(rownames(meta),function(x){unlist(strsplit(x,split=":"))[1]})
rownames(meta)=mn
id=read.table("paired.sample.txt")
id$V1=paste("S",id$V1,sep="")
id$V3=paste("S",id$V3,sep="")

sets=intersect(rownames(meta),rownames(aplicon))
subm=meta[which(rownames(meta) %in% sets),]
suba=aplicon[which(rownames(aplicon) %in% sets),]
subm=subm[rownames(suba),]

submT=subm[,which(colnames(subm) %in% id$V3)]
subaT=suba[,which(colnames(suba) %in% colnames(submT))]
subid=id[which(id$V3 %in% colnames(subaT)),]
subaC=suba[,which(colnames(suba) %in% subid$V1)]

submT=submT[,subid$V3]
subaC=subaC[,subid$V1]
subaT=subaT[,subid$V3]
pwy="PWY-5088"
m =unlist(submT[which(rownames(submT) == pwy),])
at =unlist(subaT[which(rownames(subaT) == pwy),])
ac =unlist(subaC[which(rownames(subaC) == pwy),])

library(plyr)
library(ggplot2)
library(ggrepel)
p<-cor.test(m,at,method="spearman")
data=data.frame(m,at)
p1=ggplot(data,aes(m, at)) + 
geom_point(alpha=0.8, shape=16,stat = "identity",position = "identity",color="red") + 
xlab('Fecal samples') + ylab('Tumor tissues') + 
stat_smooth(method="lm",color="black") + 
annotate(geom="text",x=8,y=5,label=paste('Spearman\'s Correlation (R^2=',round(p$estimate[[1]],3),'; p=',round(p$p.value,3),' )',sep=""),colour='black',size=3) + 
theme(axis.title.y =element_text(size=14),axis.text.y=element_text(size=12)) +  
theme(axis.title.x =element_text(size=14),axis.text.x =element_text(size=12))+
theme_classic(base_size = 13)+
ggtitle("PWY-5088") +
theme(plot.title = element_text(hjust = 0.5))

pv2<-cor.test(m,ac,method="spearman")
data2=data.frame(m,ac)
p2=ggplot(data2,aes(m, ac)) + 
geom_point(alpha=0.8, shape=16,stat = "identity",position = "identity",color="red") + 
xlab('Fecal samples') + ylab('Adjacent tissues') + 
stat_smooth(method="lm",color="black") + 
annotate(geom="text",x=8,y=5,label=paste('Spearman\'s Correlation (R^2=',round(pv2$estimate[[1]],3),'; p=',round(pv2$p.value,3),' )',sep=""),colour='black',size=3) + 
theme(axis.title.y =element_text(size=14),axis.text.y=element_text(size=12)) +  
theme(axis.title.x =element_text(size=14),axis.text.x =element_text(size=12))+
theme_classic(base_size = 13)+
ggtitle("PWY-5088") +
theme(plot.title = element_text(hjust = 0.5))


library(gridExtra)
grid.arrange(p1,p2,nrow=1)

pdf("PWY-5088.correlation.pdf")



