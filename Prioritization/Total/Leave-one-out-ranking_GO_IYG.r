load("~/prioritization/Leave-one-out/Total/Source/GO.RData")
d=G
iyg=read.table("~/prioritization/Combined_Network/Filtering_list/cregiyg.txt")
iyg=as.vector(unlist(iyg))
library(e1071)
library(randomForest)
b_pvalue=matrix(NA,length(iyg),length(d))
b_odd=matrix(NA,length(iyg),length(d))
b_Jac_Q=matrix(NA,length(iyg),length(d))
b_Jac_N=matrix(NA,length(iyg),length(d))
b_Jac=matrix(NA,length(iyg),length(d))
b_pvalue_1=matrix(NA,length(iyg),length(d))
b_odd_1=matrix(NA,length(iyg),length(d))

for (i in 1:length(d)){
  a=as.data.frame(d[i])
  a=cbind(c(1:dim(a)[1]),a)
  print(i)
  
  a_pvalue=a[order(a[,2]),]
  a_odd=a[order(-a[,3]),]
  a_Jac_Q=a[order(-a[,6]),]
  a_Jac_N=a[order(-a[,7]),]
  a_Jac=a[order(-a[,8]),]
  a_pvalue_1=a[order(a[,4]),]
  a_odd_1=a[order(-a[,5]),]
  
  for (j in 1:length(iyg)){
    	b_pvalue[j,i]=which(a_pvalue[,1]==iyg[j])
	b_odd[j,i]=which(a_odd[,1]==iyg[j])
	b_Jac_Q[j,i]=which(a_Jac_Q[,1]==iyg[j])
	b_Jac_N[j,i]=which(a_Jac_N[,1]==iyg[j])
	b_Jac[j,i]=which(a_Jac[,1]==iyg[j])
	b_pvalue_1[j,i]=which(a_pvalue_1[,1]==iyg[j])
	b_odd_1[j,i]=which(a_odd_1[,1]==iyg[j])
  }
}

com=data.frame(diag(b_pvalue),diag(b_odd),diag(b_pvalue_1),diag(b_odd_1),diag(b_Jac_Q),diag(b_Jac_N),diag(b_Jac))

write.table(com,file="~/prioritization/Leave-one-out/Total/Results/GO.txt",sep="\t",col.names=T,row.names=F)
