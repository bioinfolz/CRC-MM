library(metafor)
dataset=c("PRJEB6070","PRJEB35990","PRJNA284355","PRJNA325650","NJMU")

###EC
data=read.delim("16s.EC.matrix.txt")
ecname="EC:2.1.1.13"
j=which(rownames(data)==ecname)
ve=NULL
res=""
	for(i in dataset){
		da=read.table(paste(i,".EC.4metafor.txt",sep=""),head=T,row.names=1)
		ve=rbind(ve,da[j,])
	}
    ve$study=dataset	
my_data <- escalc(n1i = len_case, n2i = len_control, m1i = mean_case, m2i = mean_control, 
                  sd1i = sd_case, sd2i = sd_control, data = ve, measure = "SMD", 
                  append = TRUE)
ma_model_1 <- rma(yi, vi, data = my_data)
p1 <- forest(ma_model_1,slab = ve$study)

##pathway
data=read.delim("16s.path.matrix.txt")
ecname="P162-PWY"
j=which(rownames(data)==ecname)
ve=NULL
res=""
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



