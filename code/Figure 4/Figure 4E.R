data=read.delim("DSS.CpG.dmrs.annotation.xls")
loc=sapply(data$Annotation,function(x){unlist(strsplit(x,split=" "))[1]})
data$loc=loc
de=read.delim("DSS.CpG.dmrs.13.res.txt")
rownames(de)=paste(de$chr,de$start,sep=":")
rownames(data)=paste(data$Chr,data$Start,sep=":")
k=merge(de,data,by=0)
hyper=k[which(k$diff.Methy>0),]
hypo =k[which(k$diff.Methy<0),]

a=data.frame(table(hypo$loc))
b=data.frame(table(hyper$loc))
a$percent=round(a$Freq/sum(a$Freq)*100,1)
b$percent=round(b$Freq/sum(b$Freq)*100,1)

pa = ggplot(a, aes(x="", y= percent, fill= Var1))+
  geom_bar(stat="identity",width=1)+
  coord_polar(theta="y")+
  labs(x="", y="",title="Differential hyper-methylation annotation")+
  theme(axis.ticks = element_blank(), axis.text.x = element_blank())+
  scale_fill_brewer(palette="Paired",breaks = a$Var1, labels = as.character(a$Var1),direction=-1)+
  geom_text(position = position_stack(reverse =F,vjust=0.5),aes(x=1.7,label =paste(percent,"%",sep="")), show.legend = FALSE, color="black")+
  theme(plot.title = element_text(hjust = 0.5), 
        legend.title = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA), 
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        axis.line= element_blank(),
        plot.background = element_rect(fill = "transparent",colour = NA))

pb = ggplot(b, aes(x="", y= percent, fill= Var1))+
  geom_bar(stat="identity",width=1)+
  coord_polar(theta="y")+
  labs(x="", y="",title="Differential hypo-methylation annotation")+
  theme(axis.ticks = element_blank(), axis.text.x = element_blank())+
  scale_fill_brewer(palette="Paired",breaks = a$Var1, labels = as.character(a$Var1),direction=-1)+
  geom_text(position = position_stack(reverse =F,vjust=0.5),aes(x=1.7,label =paste(percent,"%",sep="")), show.legend = FALSE, color="black")+
  theme(plot.title = element_text(hjust = 0.5), 
        legend.title = element_blank(),legend.position = "bottom",
        panel.background = element_rect(fill = "transparent",colour = NA), 
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        axis.line= element_blank(),
        plot.background = element_rect(fill = "transparent",colour = NA))

				
library(grid)
grid.newpage()
pushViewport(viewport(layout=grid.layout(1,2)))
vplayout=function(x,y)
viewport(layout.pos.row=x,layout.pos.col=y)
print(pa,vp=vplayout(1,1))
print(pb,vp=vplayout(1,2))		



