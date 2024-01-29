###genus level heatmap
data=read.table("genus.raw.txt",head=T,sep="\t",row.names=1)
spdata<-read.table("sample.metadata.tsv",head=T)
data<-data[,spdata$sampleid]
data<-data+1
comp=apply(data,2,function(x){x/sum(x)})
data<-comp
res=read.delim("tissue.lefse.genus.LDA2.res.txt")
sig=res[which(res$V4>2),]
sigdf<-data[which(rownames(data) %in% sig$V1),]
library(ComplexHeatmap)
library(circlize)

log=log(sigdf)
group=c(rep("Adjacent",24),rep("Tumor",24))
Group = factor(group,levels = c("Adjacent","Tumor"))

col_fun = colorRamp2(c(-2, 0, 2), c("#2fa1dd", "white", "#f87669"))
top_annotation = HeatmapAnnotation(
cluster = anno_block(gp = gpar(fill = c("#2fa1dd", "#f87669")),
                       labels = c("Adjacent","Tumor"),
                       labels_gp = gpar(col = "white", fontsize = 12)))

m = Heatmap(t(scale(t(log))),name = " ",
            col = col_fun,
        top_annotation = top_annotation,
        show_heatmap_legend = T,
        border = F,
        show_column_names = F,
        show_row_names =T,
        column_title = NULL,
		cluster_columns = FALSE,
		column_split = Group,

		)
m


#####species
data=read.table("species.raw.txt",head=T,sep="\t")
data<-data[,spdata$sampleid]
myd=data+1
comp=apply(myd,2,function(x){x/sum(x)})
df=comp
res=read.delim("tissue.lefse.species.LDA2.res.txt")
sig=res[which(res$V4>2),]
sigdf<-df[which(rownames(df) %in% sig$V1),]
log=log(sigdf)
group=c(rep("Adjacent",24),rep("Tumor",24))
Group = factor(group,levels = c("Adjacent","Tumor"))
col_fun = colorRamp2(c(-2, 0, 2), c("#2fa1dd", "white", "#f87669"))
top_annotation = HeatmapAnnotation(
cluster = anno_block(gp = gpar(fill = c("#2fa1dd", "#f87669")),
                       labels = c("Adjacent","Tumor"),
                       labels_gp = gpar(col = "white", fontsize = 12)))

m = Heatmap(t(scale(t(log))),name = " ",
         col = col_fun,
        top_annotation = top_annotation,
        show_heatmap_legend = T,
        border = F,
        show_column_names = F,
        show_row_names =T,
        column_title = FALSE,
		cluster_columns = F,
		  column_split = Group,
		)
m
