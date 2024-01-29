gn=read.delim("dmpN_lasso_hdi.stabs.res.beta.sig.anno.txt",head=T)
gt =read.delim("dmpT_lasso_hdi.stabs.res.beta.sig.anno.txt",head=T)
gn=gn[which(gn$type=="protein-coding"),]
gt=gt[which(gt$type=="protein-coding"),]

gn$microbe=sapply(as.character(gn$taxa),function(x){unlist(strsplit(x,split="[|]"))[2]})
gt$microbe=sapply(as.character(gt$taxa),function(x){unlist(strsplit(x,split="[|]"))[2]})
gn$dire=ifelse(gn$beta >0,1,-1)
gt$dire=ifelse(gt$beta >0,1,-1)
gt2=gt[,c(19,17,5,20)]
gn2 =gn[,c(18,16,4,19)]

gt_node=c(unique(gt2$microbe),unique(gt2$gene.1))
gt_label=c(rep("microbe",length(unique(gt2$microbe))),rep("gene",length(unique(gt2$gene.1))))
gt_n =data.frame(gt_node,gt_label)

gn_node =c(unique(gn2$microbe),unique(gn2$gene.1))
gn_label =c(rep("microbe",length(unique(gn2$microbe))),rep("gene",length(unique(gn2$gene.1))))
gn_n=data.frame(gn_node,gn_label)

write.table(gt2,"tumor.genus.edge4cyto.txt",quote=F,sep="\t",row.names=F)
write.table(gt_n,"tumor.genus.node4cyto.txt",quote=F,sep="\t",row.names=F)

write.table(gn2,"normal.genus.edge4cyto.txt",quote=F,sep="\t",row.names=F)
write.table(gn_n,"normal.genus.node4cyto.txt",quote=F,sep="\t",row.names=F)


###
gn=read.delim("species.dmpN_lasso_hdi.stabs.res.beta.sig.anno.txt",head=T)
gt =read.delim("species.dmpT_lasso_hdi.stabs.res.beta.sig.anno.txt",head=T)
gn=gn[which(gn$type=="protein-coding"),]
gt=gt[which(gt$type=="protein-coding"),]
gn$microbe=gn$taxa
gt$microbe=gt$taxa
gn$dire=ifelse(gn$beta >0,1,-1)
gt$dire=ifelse(gt$beta >0,1,-1)
gt2=gt[,c(18,16,4,19)]
gn2 =gn[,c(18,16,4,19)]

gt_node=c(unique(gt2$microbe),unique(gt2$gene.1))
gt_label=c(rep("microbe",length(unique(gt2$microbe))),rep("gene",length(unique(gt2$gene.1))))
gt_n =data.frame(gt_node,gt_label)

gn_node =c(unique(gn2$microbe),unique(gn2$gene.1))
gn_label =c(rep("microbe",length(unique(gn2$microbe))),rep("gene",length(unique(gn2$gene.1))))
gn_n=data.frame(gn_node,gn_label)

write.table(gt2,"tumor.species.edge4cyto.txt",quote=F,sep="\t",row.names=F)
write.table(gt_n,"tumor.species.node4cyto.txt",quote=F,sep="\t",row.names=F)

write.table(gn2,"normal.species.edge4cyto.txt",quote=F,sep="\t",row.names=F)
write.table(gn_n,"normal.species.node4cyto.txt",quote=F,sep="\t",row.names=F)


####igraph
library(igraph）
data1 =read.table("tumor.genus.edge4cyto.txt",head=T)
data2 =read.table("normal.genus.edge4cyto.txt",head=T)
data3 =read.table("normal.species.edge4cyto.txt",head=T)
data4 =read.table("tumor.species.edge4cyto.txt",head=T)

tumor =rbind(data1,data4)
normal =rbind(data2,data3)

###
ndscale = function(x,a,b,c,d){c + (x-a)/(b-a)*(d-c)}
	library(igraph)
	library(RColorBrewer)
		m=unique(tumor$microbe)
		gg=unique(tumor$gene.1)
		col=NULL
		
		col[seq(1,length(m))]=brewer.pal(3,"Blues")[3]	  
		col[seq(length(m)+1,length(m)+length(gg))]=brewer.pal(8,"Set1")[5]
		kk=tumor[,c(1,2)]		

		weight<-tumor$r.sqr
		weight<-ndscale(weight,min(weight),max(weight),0.6,5)
		g <- graph_from_data_frame(kk, directed=FALSE)
		d<-degree(g, v=V(g), mode = c("all", "out", "in", "total"),loops = TRUE, normalized = TRUE)
		if(max(d)>min(d)){d<-ndscale(d, min(d), max(d), 6, 15)}else{d<-d/min(d)*10}
		
		V(g)$color<-col
		V(g)$size<-7
		V(g)$frame.color<-NA
		E(g)$width<-weight
		E(g)$color<-ifelse(tumor$dire >0,brewer.pal(3,"Paired")[3],brewer.pal(5,"Paired")[5]	) ##绿色positive association 红色，negative association
		V(g)$label.cex<-0.7
		V(g)$label.dist=1
		V(g)$label.color=brewer.pal(8,"RdGy")[8]
		plot(g,layout = layout.auto)

###normal

ndscale = function(x,a,b,c,d){c + (x-a)/(b-a)*(d-c)}
	library(igraph)
	library(RColorBrewer)
		m=unique(normal$microbe)
		gg=unique(normal$gene.1)
		col=NULL
		
		col[seq(1,length(m))]=brewer.pal(3,"Blues")[3]	  
		col[seq(length(m)+1,length(m)+length(gg))]=brewer.pal(8,"Set1")[5]
		kk= normal[,c(1,2)]		

		weight<-normal$r.sqr
		weight<-ndscale(weight,min(weight),max(weight),0.6,5)
		g <- graph_from_data_frame(kk, directed=FALSE)
		if(max(d)>min(d)){d<-ndscale(d, min(d), max(d), 6, 15)}else{d<-d/min(d)*10}
		
		V(g)$color<-col
		V(g)$size<-7
		V(g)$frame.color<-NA
		E(g)$width<-weight
		E(g)$color<-ifelse(normal$dire >0,brewer.pal(3,"Paired")[3],brewer.pal(5,"Paired")[5]	)
		V(g)$label.cex<-0.7
		V(g)$label.dist=1
		V(g)$label.color=brewer.pal(8,"RdGy")[8]
		plot(g,layout = layout.auto)











