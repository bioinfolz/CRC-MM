data=read.table("fecal.pathway.4lefse.res",head=F,sep="\t")
sig=data[which(data$V4>2),]
sig$path.name<-sapply(sig$V1,function(x){a=unlist(strsplit(x,split="_"));paste(a[1],a[2],sep="-")})
raw=read.delim("path.raw.txt",head=T,row.names=1)
id=sapply(rownames(raw),function(x){unlist(strsplit(x,split=":"))[1]})
raw$id=id
sig$V1<-rownames(raw[match(sig$path.name,raw$id),])
sig$V1<-sapply(sig$V1,function(x){unlist(strsplit(x,split=":"))[2]})

sig[1,1]="cholate degradation (bacteria, anaerobic)"
sig[10,1]="heme biosynthesis I (aerobic)"
sig[12,1]="methanogenesis from acetate"
write.table(sig,"meta.path.sig.txt",quote=F,sep="\t")
###lefse bar
sig<-sig[which(sig$V4>3),] ###only plot pathways with LDA>3
library("RColorBrewer")
colset <- brewer.pal(3, "Set1")
library(ggplot2)
sig$LDA=ifelse(sig$V3=="CRC",sig$V4,-sig$V4)
sig<-sig[order(sig$LDA),]
sig$V1=factor(sig$V1,levels=(as.character(sig$V1)))
library(ggplot2)
ggplot(sig,aes(x=V1,y=LDA,fill=as.character(sig$V3))) +
geom_bar(stat="identity",position="dodge")+
scale_fill_manual(values= colset[c(1,2)] )+
coord_flip()+xlab("")+theme_classic()+
theme(legend.title=element_blank())+
geom_hline(yintercept = c(-4,-3,-2,-1),linetype =2)+
geom_hline(yintercept = c(1,2,3,4),linetype =2)

