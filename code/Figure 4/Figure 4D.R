library(stringr)
info= read.delim("TCGA-COAD.GDC_phenotype.tsv",head=T,sep="\t")
pd=info
pd$sampleID = str_sub(pd$submitter_id.samples,1,15)
pd$patient = str_sub(pd$sampleID,1,12)
pd$group_list = ifelse(as.numeric(str_sub(pd$sampleID,14,15))<10,"Tumor","Normal")
table(pd$group_list)

tp = pd$patient[pd$group_list =="Normal"]
np = pd$patient[pd$group_list =="Tumor"]
okp = intersect(tp,np)

pd <- pd[pd$patient %in% okp,]
label<-str_sub(pd$submitter_id.samples,16,16)
pd<-pd[which(label=="A"),]
rownames(pd) <- pd$sampleID
pd<-pd[,c(121,122,123)]
pd1=pd

library(stringr)
info= read.delim("TCGA-READ.GDC_phenotype.tsv",head=T,sep="\t")
pd=info
pd$sampleID = str_sub(pd$submitter_id.samples,1,15)
pd$patient = str_sub(pd$sampleID,1,12)
pd$group_list = ifelse(as.numeric(str_sub(pd$sampleID,14,15))<10,"Tumor","Normal")
table(pd$group_list)

tp = pd$patient[pd$group_list =="Normal"]
np = pd$patient[pd$group_list =="Tumor"]
okp = intersect(tp,np)

pd <- pd[pd$patient %in% okp,]
label<-str_sub(pd$submitter_id.samples,16,16)
pd<-pd[which(label=="A"),]
rownames(pd) <- pd$sampleID
pd<-pd[,c(119,120,121)]

df=rbind(pd1,pd)
write.table(df,"CRC.pair.group.list.txt",quote=F,sep="\t")


library(data.table)
b = data.table::fread("methylation/HumanMethylation450.gz", data.table = F)
a <-b
rownames(a)=a$sample

d = a[,colnames(a) %in% df$sampleID]
patient = str_sub(colnames(d),1,12)
group_list = ifelse(as.numeric(str_sub(colnames(d),14,15))<10,"Tumor","Normal")

tp = patient[group_list =="Normal"]
np = patient[group_list =="Tumor"]
okp = intersect(tp,np);length(okp)
d = d[,patient %in% okp];dim(d)

df = df[match(colnames(d),df$sampleID),]
identical(str_sub(colnames(d),1,15),df$sampleID)

G1=d[,colnames(d) %in% df[which(df$group_list=="Tumor"),1]]
G2=d[,colnames(d) %in% df[which(df$group_list=="Normal"),1]]
dim(G1);dim(G2)
d=cbind(G1,G2)
write.table(d,"CRC.450K.CC.matrix.txt",quote=F,sep="\t")
save(d,file="CRC.CC.matrix.Rdata")


load("CRC.CC.matrix.Rdata")
beta=as.matrix(d)
beta=impute.knn(beta) 
sum(is.na(beta))
beta=beta$data
beta=beta+0.00001
save(beta,file="CRC.paired.450K.beta.Rdata")

library(ChAMP)
load("CRC.paired.450K.beta.Rdata")
sampleID=colnames(beta)
pd=read.table("CRC.pair.group.list.txt",head=T)
pd=pd[which(pd$sampleID %in% colnames(beta)),]
beta=beta[,pd$sampleID] 
myLoad=champ.filter(beta = beta,pd=pd) #这一步已经自动完成了过滤
dim(myLoad$beta)
save(myLoad,file = 'step1_myLoad.Rdata')

norm_file = "./step2_champ_myNorm.Rdata"
myNorm <- champ.norm(beta=myLoad$beta,arraytype="450K",cores=2)
save(myNorm,file = norm_file)


load("step2_champ_myNorm.Rdata")
group_list <- pd$group_list

group_list<-rep(c("Tumor","adjacent"),45)
myDMP <- champ.DMP(beta = myNorm,pheno=group_list)
myDMR <- champ.DMR(beta = myNorm,pheno=group_list,method="Bumphunter")
save(myDMR,file="CRC.paired.90.myDMR.res.Rdata")
save(myDMP,file="CRC.paired.90.myDMP.res.Rdata")

myGSEA<-champ.GSEA(beta=myNorm,
        DMP=myDMP[[1]],
		DMR=myDMR,
		CpGlist=NULL,
		Genelist=NULL,
		pheno=group_list,
		method="gometh",
		arraytype="450K",
		Rplot=TRUE,
		adjPval=0.01)
save(myDMR,myDMP,myGSEA,file="CRC.paired.myDMR.myDMP.myGSEA.Rdata")

res<-myGSEA$DMR
res2<-myGSEA$DMP
BP=res[which(res$ONTOLOGY=="BP"),]
MF=res[which(res$ONTOLOGY=="MF"),]


df=data.frame(myDMR)
DMR.GUI(DMR=myDMR,beta=myNorm,pheno=group_list)


library(missMethyl)
load("CRC.paired.myDMR.myDMP.myGSEA.Rdata")
library(GenomicRanges)
singleGRange <- GRanges(as.data.frame(GRangesList(myDMR)))
singleGRange <- GRanges(as.data.frame(GRangesList(myDMP)))


GO<-(
singleGRange,
all.cpg = NULL,
collection = c("GO", "KEGG"),
array.type = c("450K", "EPIC"),
plot.bias = FALSE,
prior.prob = TRUE,
anno = NULL,
equiv.cpg = TRUE,
fract.counts = TRUE,
genomic.features = c("ALL", "TSS200", "TSS1500", "Body", "1stExon", "3'UTR", "5'UTR","ExonBnd"),
sig.genes = TRUE
)

sig =GO[which(GO$FDR<0.01),]
sig=sig[order(sig$FDR),]
BP= sig[which(sig$ONTOLOGY=="BP"),]
MF =sig[which(sig$ONTOLOGY=="MF"),]
write.table(BP,"missMethyl.DMR.BP.sig.txt",quote=F,sep="\t")
write.table(MF,"missMethyl.DMR.MF.sig.txt",quote=F,sep="\t")


GO2<-goregion(
singleGRange,
all.cpg = NULL,
collection = c("GO", "KEGG"),
array.type = c("450K", "EPIC"),
plot.bias = FALSE,
prior.prob = TRUE,
anno = NULL,
equiv.cpg = TRUE,
fract.counts = TRUE,
genomic.features = c("TSS200", "TSS1500","5'UTR"),
sig.genes = TRUE
)

sig =GO2[which(GO2$FDR<0.01),]
sig=sig[order(sig$FDR),]
BP= sig[which(sig$ONTOLOGY=="BP"),]
MF =sig[which(sig$ONTOLOGY=="MF"),]
write.table(BP,"missMethyl.DMR.promoter.BP.sig.txt",quote=F,sep="\t")
write.table(MF,"missMethyl.DMR.promoter.MF.sig.txt",quote=F,sep="\t")

sig3 =GO2[which(GO2$FDR<0.0001),]
sig3=sig3[order(sig3$FDR),]
sig3=sig3[which(sig3$ONTOLOGY %in% c("MF","BP")),]

###get DMP loci for each sample
load("step2_champ_myNorm.Rdata")
str(myNorm)

load("CRC.paired.90.myDMP.res.Rdata")
sig=data.frame(myDMP)

sigdata=myNorm[which(rownames(myNorm) %in% rownames(sig)),]
write.table(sigdata,"champ.DMP.matrix.txt",quote=F,sep="\t")

promoter region:
sigp=sig[which(sig$Tumor_to_adjacent.feature %in% c("TSS1500","TSS200","5'UTR")),]
sigdata=myNorm[which(rownames(myNorm) %in% rownames(sigp)),]
write.table(sigdata,"champ.DMP.promoter.matrix.txt",quote=F,sep="\t")

###
matrix=read.table("champ.DMP.promoter.matrix.txt")
group_list<-rep(c("Tumor","adjacent"),45)
names(group_list)=colnames(matrix)

mat=matrix
mat[is.na(mat)]<-0

library(ggfortify)

pca<-prcomp(t(mat)) 
group_list<-rep(c("Tumor","adjacent"),45)

scores <- data.frame(group_list, pca$x[,1:2])

summ = summary(pca)
xlab = paste0("PC1(",round(summ$importance[2,1]*100,2),"%)")
ylab = paste0("PC2(",round(summ$importance[2,2]*100,2),"%)")


ggplot(data = scores,aes(x=PC1,y= PC2,color=group_list))+
stat_ellipse(type="norm",geom="polygon",alpha=0.1,color=NA,aes(fill=group_list))+
geom_point()+
scale_fill_manual(values=c("blue","red"))+
scale_color_manual(values=c("blue","red"))+
theme_classic()+
labs(x=xlab,y=ylab,color="")+guides(fill=F)



