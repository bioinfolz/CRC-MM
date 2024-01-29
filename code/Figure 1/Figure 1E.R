########meta analysis
samp =read.delim("meta.txt",head=T,row.names=2)
EC=read.delim("sample.count.txt",head=T,check.names=FALSE,row.names=1)
EC=EC[,-1]
comp=apply(EC,2,function(x){x/sum(x)})
asx=asin(sqrt(comp)) ##arcsine-square-root transformed

dataset=unique(as.character(samp$dataset_name))
library(metafor)
j=283
ecname=rownames(EC)[j]
ve=NULL
	for(i in dataset){
		da=read.table(paste(i,".pathway.4metafor.txt",sep=""),head=T,row.names=1)
		ve=rbind(ve,da[j,])
	}
    ve$study=dataset

my_data <- escalc(n1i = len_case, n2i = len_control, m1i = mean_case, m2i = mean_control, 
                  sd1i = sd_case, sd2i = sd_control, data = ve, measure = "SMD", 
                  append = TRUE)
ma_model_1 <- rma(yi, vi, data = my_data)
p1 <- forest(ma_model_1,slab = ve$study)
 
 
 
 
 
 
 
 
 