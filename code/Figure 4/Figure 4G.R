data=read.delim("DMR.TSS.res.pc.genelist.David.direct.txt")
BP=data[which(data$Category=="GOTERM_BP_DIRECT"),]
MF=data[which(data$Category=="GOTERM_MF_DIRECT"),]

data=data[order(data$PValue),]
sub =data[seq(1,30),c(1,2,5)]
GObarplot(sub)
