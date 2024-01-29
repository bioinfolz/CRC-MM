library(methylKit)
load("meth.13.RData")
mat=percMethylation(meth,chunk.size = 1e+06)
cs=apply(mat,2,mean)
control=cs[seq(1,26,by=2)]
tumor=cs[seq(2,26,by=2)]
df=data.frame(control,tumor)
write.table(df,"P13.percMethylation.stat.txt",quote=F,sep="\t")
##plot
library(ggplot2)
library("RColorBrewer")
colset <- brewer.pal(6, "Set1")
library(ggsignif)
my_comparisons=list(c("Normal","Tumor"))

data=read.table("P13.percMethylation.stat.txt",head=T)
Percent=c(data$control,data$tumor)
Label=rep(c("Normal","Tumor"),each=13)
Sample=rep(seq(1,13),2)
df=data.frame(Percent,Label,Sample)

pdf("percentage.paired.P13.pdf",width=4,height=8)