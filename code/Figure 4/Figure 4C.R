data=read.delim("CRC.450K.case.control.beta.msd.txt")
library(stringr)
label<-str_sub(rownames(data),1,12)

type<-str_sub(rownames(data),14,15)
data$group=ifelse(type=="01","Tumor","Adjacent")

library(ggplot2)
library(ggbeeswarm)
library(RColorBrewer)
colset <- brewer.pal(3,"Set1")
ggplot(data,aes(x=type,y=g))+
geom_beeswarm(aes(color=type),cex=2)+  scale_color_manual(values= colset[c(2,1)])+
theme_classic()+
theme(legend.position = c("top"))+ 
xlab("")

library(ggsignif)
my_comparisons=list(c("Adjacent","Tumor"))

data$patient=label
a=label[duplicated(label)]
sub=data[which(data$patient %in% a),]
my_comparisons=list(c("Adjacent","Tumor"))
ggplot(sub,  aes(x=group,y=g))+
geom_violin(aes(fill=group))+
geom_boxplot(width=0.1)+
geom_signif(color="black",comparisons=my_comparisons,
map_signif_level=T,test=wilcox.test)+
theme_classic()+
scale_fill_manual(values= colset[c(2,1)] )





