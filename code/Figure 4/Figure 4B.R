
mat=read.table("DML.CpG.TSS.deta20.NAfil.matrix.txt")

pca<-prcomp(t(mat)) 
group_list<-rep(c("Adjacent","Tumor"),13)
scores <- data.frame(group_list, pca$x[,1:2])

summ = summary(pca)
xlab = paste0("PC1(",round(summ$importance[2,1]*100,2),"%)")
ylab = paste0("PC2(",round(summ$importance[2,2]*100,2),"%)")

pdf("pca.pdf")
ggplot(data = scores,aes(x=PC1,y= PC2,color=group_list))+
stat_ellipse(type="norm",geom="polygon",alpha=0.1,color=NA,aes(fill=group_list))+
geom_point()+
scale_fill_manual(values=c("blue","red"))+
scale_color_manual(values=c("blue","red"))+
theme_classic()+
labs(x=xlab,y=ylab,color="")+guides(fill=F)
dev.off()



