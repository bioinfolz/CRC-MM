source("F://script/functions/ggplot-自定义background.R")
library("RColorBrewer")
colset <- brewer.pal(5, "Set1")
library("ggplot2")
library("vegan")
library("GUniFrac")
library("ggpubr")
library(ape)

Otu_tab <- read.table('table.from_biom.txt',row.names=1,head=T,sep="\t",check.names=FALSE)
Tree   <- read.tree('unrooted.tree.nwk')
Otu_tab<-Otu_tab[,-which(colnames(Otu_tab) %in% c("S56"))]  ##blank control
comp=apply(Otu_tab,2,function(x){x/sum(x)})
ind=apply(comp,1,function(x){length(x[which(x>0.0001)])})
Otu_tab<-Otu_tab[which(ind>4),]
Otu_tab1 <-as.data.frame(t(Otu_tab))
Otu_tab_rff<-Rarefy(Otu_tab1)$otu.tab.rff

# Variance adjusted weighted UniFrac
Tree   <- read.tree('rooted.tree.nwk')
unifracs <- GUniFrac(Otu_tab_rff,Tree,alpha=c(0,0.5,1))
du <- unifracs$unifracs[,,"d_UW"]

a<-adonis(du~sample$type,permutations = 999,method="bray")
0.03934  0.008 
 statistic R2: 0.0667
      Significance: 0.001 
	  
PCOA <- pcoa(du)
result <-PCOA$values[,"Relative_eig"]
pro1 = as.numeric(sprintf("%.3f",result[1]))*100
pro2 = as.numeric(sprintf("%.3f",result[2]))*100
pc = PCOA$vectors[,c(1,2)]
pc=data.frame(pc,sample$type)
xlab=paste("PCOA1(",pro1,"%)",sep="") 
ylab=paste("PCOA2(",pro2,"%)",sep="")
p2 =ggscatter(pc,x="Axis.1",y="Axis.2",color="sample.type",ellipse=TRUE,size=4,palette=colset[c(2,1)])+
labs(x=xlab,y=ylab,title="UniFrac PCoA")+theme_zg()+
geom_text(aes(label=rownames(points)),size=2,vjust=-1) +
 annotate("text",x=0.004,y=0.3,parse=TRUE,size=4,label="'R2='* 0.0667",family="serif",fontface="italic",colour="black")+
 annotate("text",x=0,y=0.25,parse=TRUE,size=4,label="'p='* 0.001",family="serif",fontface="italic",colour="black")
