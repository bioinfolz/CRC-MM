##TCGA tumor
data=read.delim("tumorg.list.DAVID.txt",head=T,sep="\t")
sig=data[which(data$PValue<0.05),]
sig$log=-log10(sig$PValue)
sig=sig[seq(1,20),]
sig$Term=factor(sig$Term,levels=rev(sig$Term))
pp = ggplot(sig,aes(log,Term))
pp + geom_point()
pp + geom_point(aes(size=(Fold.Enrichment),color=Category))+ylab("-log10(Pvalue)")+
xlab("-log10(Pvalue)")+theme_bw()

##TCGA normal
data=read.delim("normalg.list.DAVID.txt",head=T,sep="\t")
sig=data[which(data$PValue<0.05),]
sig$log=-log10(sig$PValue)
sig$Term=sapply(sig$Term,function(x){unlist(strsplit(x,split="~"))[2]})
sig=sig[seq(1,20),]
sig$Term=factor(sig$Term,levels=rev(sig$Term))
pp = ggplot(sig,aes(log,Term))
pp + geom_point()
pp + geom_point(aes(size=(Fold.Enrichment),color=Category))+ylab("-log10(Pvalue)")+
xlab("-log10(Pvalue)")+theme_bw()

