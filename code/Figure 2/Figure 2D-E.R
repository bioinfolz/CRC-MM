###lefse bar

sig=read.table("tissue.path.LDA2.res.txt",head=T)
anno=read.delim("path.anno.txt",head=T)
sig$id=sapply(sig$V1,function(x){gsub("_", "-", x, fixed = TRUE)})
sig$des<-anno[match(sig$id,anno$pathway),2]

library("RColorBrewer")
colset <- brewer.pal(3, "Set1")
sig$LDA=ifelse(sig$V3=="T",sig$V4,-sig$V4)  ##T for tumor, C for control
sig<-sig[order(sig$LDA),]
sig$des=factor(sig$des,levels=(as.character(sig$des)))
library(ggplot2)
p1=ggplot(sig,aes(x=des,y=LDA,fill=as.character(sig$V3))) +
geom_bar(stat="identity",position="dodge")+
scale_fill_manual(values= colset[c(2,1)] )+
coord_flip()+xlab("")+theme_classic()+
theme(legend.title=element_blank())+
geom_hline(yintercept = c(-3,-2,-1),linetype =2)+
geom_hline(yintercept = c(3,2,1),linetype =2)
write.table(sig,"tissue.path.LDA2.res.anno.txt",row.names=F,quote=F,sep="\t")

sig=read.table("tissue.EC.LDA2.res.txt",head=T)
anno=read.delim("EC.anno.txt",head=T)
sig$id=sapply(sig$V1,function(x){gsub("_", ".", x, fixed = TRUE)})
sig$des<-anno[match(sig$id,anno$function.),2]

library("RColorBrewer")
colset <- brewer.pal(3, "Set1")
sig$LDA=ifelse(sig$V3=="T",sig$V4,-sig$V4)
sig<-sig[order(sig$LDA),]
sig<-sig[which(abs(sig$LDA)>2.5),] ##plot only LDA>2.5
sig$des=factor(sig$des,levels=(as.character(sig$des)))
library(ggplot2)
p2=ggplot(sig,aes(x=des,y=LDA,fill=as.character(sig$V3))) +
geom_bar(stat="identity",position="dodge")+
scale_fill_manual(values= colset[c(2,1)] )+
coord_flip()+xlab("")+theme_classic()+
theme(legend.title=element_blank())+
geom_hline(yintercept = c(-3,-2,-1),linetype =2)+
geom_hline(yintercept = c(3,2,1),linetype =2)

