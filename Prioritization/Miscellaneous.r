

# Forgetting to calculate direct edges and shared genes with cregiyg
# then have to compute those separately and combined with other data
# This scripts combines the previous data with new one
load("text.RData")
list=read.table("~/prioritization/New_IYG_Results/order.txt")
list=as.vector(unlist(list))
a=lapply(list,read.table,sep="\t",header=T)
for(i in 1:147){
b=data.frame(w[i])
c=data.frame(a[i])
w[i]=list(data.frame(b[,1],c,b[,2:8]))
}

for(i in 1:147){
b=data.frame(w[i])
colnames(b)=c("alldeg",colnames(b)[-1])
w[i]=list(b)
}
head(data.frame(w[1]))
save(w,file="text.RData")







list=read.table("~/prioritization/New_IYG_Results/order.txt")
list=as.vector(unlist(list))
a=lapply(list,read.table,sep="\t",header=T)
q=matrix(NA,1,dim(data.frame(a[1]))[2])
colnames(q)=colnames(data.frame(a[1]))

for (i in 1:147){
b=data.frame(a[i])
q=rbind(q,b[i,])
}
q=q[-1,]
head(q)
apply(q,2,quantile)
write.table(q,file="PCC_Results.txt",row.names=F,sep="\t",quote=F)







## For all supervised and ensemble approaches - modification needed for removing some classifiers
# Mapping back from Networks codes to AT
a=read.table("~/prioritization/Ranking/Total.txt",header=T)
a=a[,-36]
r=read.table("~/prioritization/Combined_Network/Ref_new.txt",sep="\t")
q=matrix(NA,dim(a)[1],c(dim(a)[2]-1))
n=data.frame(NA)
for (j in 2:dim(a)[2]){
	print(j)
	if (j%in%c(11,12,13,14,18,19,20,24,25,26,30,31,32,51,52,53,54,55,56,57,58,59,60)) b=a[order(a[,j]),1]
	else b=a[order(-a[,j]),1]
	for (i in 1:200){
		q[i,j-1]=as.vector(unlist(r[which(r[,2]==b[i]),1]))
	}
	n=data.frame(n,b)
}
n=n[,-1]

# Mapping back from Networks AT to gene's names
r=read.csv("~/prioritization/Combined_Network/TAIR/REF_New20130130.txt",sep="\t")
d=matrix(NA,dim(q)[1],dim(q)[2])
for (j in 1:dim(q)[2]){
	for (i in 1:dim(q)[1]){
		if(length(as.vector(unlist(r[which(r[,1]==q[i,j]),2])))==0) d[i,j]=q[i,j]
		else d[i,j]=as.vector(unlist(r[which(r[,1]==q[i,j]),2]))
	}
}


y=read.table("prioritization/Combined_Network/New/Filtering_list/y.txt")
for (j in 1:dim(q)[2]){
	for (i in 1:200){
		if(length(as.vector(unlist(r[which(r[,1]==q[i,j]),2])))==0) d[i,j]=q[i,j]
		else d[i,j]=as.vector(unlist(r[which(r[,1]==q[i,j]),2]))
	}
}

y=read.table("prioritization/Combined_Network/New/Filtering_list/y.txt")
y1=matrix(NA,dim(a)[1],dim(a)[2]-1)
for (j in 2:3){
	for (i in 1:dim(a)){
		y1[i,j]=y[a[i,j],]
	}
}



colnames(q)=colnames(d)=colnames(n)=colnames(a)[-1]
write.table(n,file="~/prioritization/Ranking/Net_code.txt",row.names=F,sep="\t",quote=F)
write.table(q,file="~/prioritization/Ranking/AT.txt",row.names=F,sep="\t",quote=F)
write.table(d,file="~/prioritization/Ranking/Name.txt",row.names=F,sep="\t",quote=F)

freq=sort(table(d[,35:dim(d)[2]]))
write.table(freq,file="~/prioritization/Ranking/High_freq.txt",sep="\t",quote=F,col.names=F)


for (j in 1:dim(c)[2]){
	for (i in 1:100){
	#if length(which(x[i,1]==c[,j]))==0 
	ra[i,j]=which(x[i,1]==c[,j])
	}
}




boxplot(t$RF,t$RF.SVM.LDA,t$RF.LDA,t$RF.NaiveBayes.LDA,t$RF.SVM.Glmnet,t$RF.SVM.LDA.NaiveBayes,t$RF.SVM.LDA.NaiveBayes.Glmnet,t$RF.LDA.Glmnet,t$RF.Glmnet,t$LDA.Glmnet,t$ens1,t$ens2,t$ens3,t$ens4,t$ens5,t$ens6,t$ens7,t$ens8,t$ens9,t$ens10,names=c("R","R.S.L","R.L","R.N.L","R.S.G","R.S.L.N","R.S.L.N.G","R.L.G","R.G","L.G","ens1","ens2","ens3","ens4","ens5","ens6","ens7","ens8","ens9","ens10"),ylab="GR genes ranking",las=2)
names=c("RF","RF.SVM.LDA","RF.LDA","RF.NaiveBayes.LDA","RF.SVM.Glmnet","RF.SVM.LDA.NaiveBayes","RF.SVM.LDA.NaiveBayes.Glmnet","RF.LDA.Glmnet","RF.Glmnet","LDA.Glmnet","ens1","ens2","ens3","ens4","ens5","ens6","ens7","ens8","ens9","ens10"),ylab="GR genes ranking"


























