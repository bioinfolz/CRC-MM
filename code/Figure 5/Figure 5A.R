###pheatmap
gn<-read.delim("control.DSS.global.cor.res.txt")
gt<-read.delim("tumor.DSS.global.cor.res.txt")
sn<-read.delim("species.control.DSS.global.cor.res.txt")
st<-read.delim("species.tumor.DSS.global.cor.res.txt")

gmerge<-merge(gn,gt,by=0)
gsig<-gmerge[which(gmerge$V1.x<0.05 | gmerge$V1.y <0.05),]
gp<-cbind(gsig[,2],gsig[,5])
gv<-cbind(gsig[,3],gsig[,6])
id<-sapply(as.character(gsig$Row.names),function(x){unlist(strsplit(x,split="[|]"))[2]})
rownames(gp)<-id
rownames(gv)<-id
colnames(gp)<-c("AN","Tumor")
colnames(gv)<-c("AN","Tumor")

##species
smerge<-merge(sn,st,by=0)
ssig <-smerge[which(smerge$V1.x<0.05 | smerge$V1.y <0.05),]
sp<- cbind(ssig[,2],ssig[,5])
sv<- cbind(ssig[,3],ssig[,6])
id<-sapply(as.character(ssig$Row.names),function(x){unlist(strsplit(x,split="[|]"))[2]})

rownames(sp)=id
rownames(sv)=id
colnames(sp)=c("AN","Tumor")
colnames(sv)=c("AN","Tumor")

library(gplots)

Mycol <-colorpanel(50,"forestgreen","white",high="#B03060") 
noteTrans <- function(matrix,C2,C3  ){
 P_list <- as.vector(matrix)
 i <- 1
 note <- c()	
 	while (i <= length(P_list)){
 		note[i] = ""
 		if (is.na(P_list[i]) == TRUE){
 		P_list[i] = 0	
 		}

 			if (P_list[i] <= C2){
 				note[i] = "*"
 			}
 			if (P_list[i] <= C3){
 				note[i] = "**"
 			}
 		i = i + 1	
 	}
 noteM <- matrix(note,byrow=F,nrow=dim(matrix)[1])
 return(noteM)	
 }



MyNote <- noteTrans(gp ,0.05,0.01)
p1=heatmap.2(gv, trace="none",density="none",col=Mycol,  offsetCol = -0.8,offsetRow = -0.4,hclustfun = function(c)hclust(c,method="ward.D2"),
cexRow=1, cexCol=1,  srtCol=45,  cellnote= MyNote, notecex=1,  notecol="black", margins=c(10,10))


MyNote <- noteTrans(sp ,0.05,0.01)
p2=heatmap.2(sv, trace="none",density="none",col=Mycol,  offsetCol = -0.8,offsetRow = -0.4,hclustfun = function(c)hclust(c,method="ward.D2"),
cexRow=1, cexCol=1,  srtCol=45,cellnote= MyNote, notecex=1,  notecol="black", margins=c(10,10))

library(ggpubr)
ggarrange(p2,p1, ncol=2, nrow=1, common.legend = TRUE, legend="bottom") 



library(grid)
grid.newpage()
pushViewport(viewport(layout=grid.layout(1,2)))
vplayout=function(x,y)
viewport(layout.pos.row=x,layout.pos.col=y)
print(p1,vp=vplayout(1,1))
print(p2,vp=vplayout(1,2))

library(gridExtra)
grid.arrange(p1, p2, ncol=2, nrow=1)






