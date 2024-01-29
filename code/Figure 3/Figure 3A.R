setwd("E://1-CRC/correlation/overall-correlation")
ad=read.delim("meta.adjacent.function.paired.correlation.res.txt")
tumor=read.delim("meta.tumor.function.paired.correlation.res.txt")
ad2=ad[,-1]
colnames(ad2)=c("ID","R.ad","P.ad","D.ad")
t2= tumor[,-1]
colnames(t2)=c("ID","R.tumor","P.tumor","D.tumor")


mg=full_join(ad2,t2,by= "ID")
mg$anno=ifelse(is.na(mg$D.ad),mg$D.tumor,mg$D.ad)
mg2 = mg[,c(1,3,4,7,8,9)]

rownames(mg)=mg$anno
mg2 = mg[,c(2,5)]
mg2[is.na(mg2)]<- 0
pheatmap(mg2,treeheight_row=0,treeheight_col=0,gaps_row=c(13,15) )
