tumor =read.table("Tumor.location.4chiq.txt",head=T) 
normal =read.table("Normal.location.4chiq.txt",head=T) 

td=data.frame(tumor[,c(2,4)],485577,3671)
td$Freq=td$Freq-1
out=NULL
for(i in seq(1,6)){
x=unlist(td[i,])
p =phyper(x[1],x[2],x[3]-x[2],x[4], lower.tail=T)
out=c(out,p)
}
res=data.frame(tumor,out)
res$enrich=res$Rcor/res$Rbg


td=data.frame(normal[,c(2,4)],485577,1484)
td$Freq=td$Freq-1
out=NULL
for(i in seq(1,6)){
x=unlist(td[i,])
p =phyper(x[1],x[2],x[3]-x[2],x[4], lower.tail=T)
out=c(out,p)
}
res2=data.frame(normal,out)
res2$enrich=res2$Rcor/res2$Rbg


res2=res2[order(res2$enrich),]
res2$k=factor(res2$k,levels=rev(res2$k))
res2$log=-log10(res2$out)
p2=ggplot(res2,aes(x=enrich,y=k,fill=log))+
geom_bar(stat="identity",width=0.7,color="navy")+ 
theme_classic()+scale_x_reverse()+
theme(axis.ticks.y = element_blank()) +
geom_vline (xintercept = c(1),linetype ="dotted",color="red")

rownames(res2)=res2$k
rownames(res)=res$k
res=res[rownames(res2),]

res$k=factor(res$k,levels=rev(res$k))
res$log=-log10(res$out)
p=ggplot(res,aes(x=enrich,y=k,fill=log))+
geom_bar(stat="identity",width=0.7,color="navy")+ 
theme_classic()+
theme(axis.ticks.y = element_blank()) +
geom_vline (xintercept = c(1),linetype ="dotted",color="red")

library(ggpubr)
ggarrange(p2,p, ncol=2, nrow=1, common.legend = TRUE, legend="bottom") 











