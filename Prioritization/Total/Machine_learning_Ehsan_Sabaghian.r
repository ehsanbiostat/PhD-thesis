
# Reading data
load("~/prioritization/Leave-one-out/Total/Source/DATA(2).RData")  # Whole Network data sets
cregiyg=read.table("~/prioritization/Combined_Network/Filtering_list/cregiyg.txt") # Known genes list (IYG & REG)
load("~/prioritization/Leave-one-out/Total/Source/GO.RData")  # GO analysis

# Upload packages
library(e1071)
library(randomForest)
library(gbm)
library(Matrix)
library(igraph)
library(MASS)

# Creating response vector
cregiyg=unlist(cregiyg)
y=sparseMatrix(c(cregiyg,27290),rep(1,147),x=1)
y[27290,]=0
y=as.matrix(y)


# Creating ranking matrix
b_bet=matrix(NA,length(cregiyg),length(d))
b_deg=matrix(NA,length(cregiyg),length(d))
b_IYG=matrix(NA,length(cregiyg),length(d))
b_IYG_S=matrix(NA,length(cregiyg),length(d))
b_REG=matrix(NA,length(cregiyg),length(d))
b_REG_S=matrix(NA,length(cregiyg),length(d))
b_clos=matrix(NA,length(cregiyg),length(d))
b_klein=matrix(NA,length(cregiyg),length(d))
b_sim.jacard=matrix(NA,length(cregiyg),length(d))
b_sim.dic=matrix(NA,length(cregiyg),length(d))
b_sim.we=matrix(NA,length(cregiyg),length(d))
b_short=matrix(NA,length(cregiyg),length(d))
b_pvalue_Top=matrix(NA,length(cregiyg),length(d))
b_pvalue_first=matrix(NA,length(cregiyg),length(d))
b_Jaccard_Q=matrix(NA,length(cregiyg),length(d))
b_Jaccard_N=matrix(NA,length(cregiyg),length(d))
b_Jaccard=matrix(NA,length(cregiyg),length(d))
b_NB=matrix(NA,length(cregiyg),length(d))
b_RF=matrix(NA,length(cregiyg),length(d))
b_gbm=matrix(NA,length(cregiyg),length(d))
b_log=matrix(NA,length(cregiyg),length(d))
b_svm=matrix(NA,length(cregiyg),length(d))
b_lda=matrix(NA,length(cregiyg),length(d))
b_ens=matrix(NA,length(cregiyg),length(d))



# Reading each left-out dataset
for (i in 1:length(cregiyg)){
  print(i)
  a=as.data.frame(d[i])
  b=as.data.frame(G[i])
  a=cbind(c(1:dim(a)[1]),a,b[,-c(2,4)],y)
  a1=a[-cregiyg[i],]
    
  # Datasets
  new=test=data.frame(a[cregiyg[i],-c(1,dim(a)[2])])
  x=data.frame(a1[,-c(1,dim(a1)[2])])
  Y=data.frame(a1[,dim(a)[2]])
  
  # Naive Bayes
  NB=naiveBayes(x,unlist(Y))
  p=predict(NB,new,type=c("raw"))[2]
  fit=predict(NB,x,type=c("raw"))[,2]
  a=cbind(a,c(fit[1:c(cregiyg[i]-1)],p,fit[cregiyg[i]:length(fit)]))
  colnames(a)[dim(a)[2]]=c("NaiveBayes")
    
  # Random Forest
  rf=randomForest(x,as.factor(unlist(Y)))
  fit=rf$votes[,2]
  p=predict(rf,new,type="vote")[2]
  a=cbind(a,c(fit[1:c(cregiyg[i]-1)],p,fit[cregiyg[i]:length(fit)]))
  colnames(a)[dim(a)[2]]=c("RF")
  
  # GBM
  gb=gbm.fit(x,unlist(Y),distribution="adaboost",n.trees=100)
  p=predict(gb,new,n.trees=100)
  fit=gb$fit
  a=cbind(a,c(fit[1:c(cregiyg[i]-1)],p,fit[cregiyg[i]:length(fit)]))
  colnames(a)[dim(a)[2]]=c("GBM")
  
     
  # Logistic Regression
  data=data.frame(Y=a1[,14],Degree=a1[,2],IYG=a1[,3],IYG_shared=a1[,4],REG=a1[,5],REG_shared=a1[,6],Betwee=a1[,7],clos=a1[,8]
  ,klein=a1[,9],sim.jacard=a1[,10],sim.dic=a1[,11],sim.we=a1[,12],short=a1[,13],pvalue_TopTen=a1[,14],pvalue_first=a1[,15],Jaccard_Q=a1[,16],Jaccard_N=a1[,17],Jaccard=a1[,18])
  model=glm(Y~Degree+IYG+IYG_shared+REG+REG_shared+Betwee+clos+klein+sim.jacard+sim.dic+sim.we+short
  +pvalue_TopTen+pvalue_first+Jaccard_Q+Jaccard_N+Jaccard,data=data,family="binomial")
  newglm=a[cregiyg[i],c(2:18)]
  p=predict(model,newglm,type=c("response"))
  fit=predict(model,type=c("response"))
  a=cbind(a,c(fit[1:c(cregiyg[i]-1)],p,fit[cregiyg[i]:length(fit)]))
  colnames(a)[dim(a)[2]]=c("Logistic")
  
  # SVM
  sv=svm(x,Y,probability=TRUE)
  fit=sv$fitted
  p=predict(sv,new,type="response")
  a=cbind(a,c(fit[1:c(cregiyg[i]-1)],p,fit[cregiyg[i]:length(fit)]))
  colnames(a)[dim(a)[2]]=c("SVM")
  
  # Linear Discriminante Analysis (Fisher)
  LDA=lda(x[,-7],unlist(Y))
  fit=predict(LDA,method=c("predictive"))$posterior[,2]
  p=predict(LDA,new[,-7],method=c("predictive"))$posterior[,2]
  a=cbind(a,c(fit[1:c(cregiyg[i]-1)],p,fit[cregiyg[i]:length(fit)]))
  colnames(a)[dim(a)[2]]=c("LDA")
  
  # Ensemble
  Ens=a$RF+a$Logistic+a$LDA
  a=cbind(a,Ens)
  
  # Ranking based on different criteria
  a_deg=a[order(-a[,2]),]
  a_IYG=a[order(-a[,3]),]
  a_REG=a[order(-a[,4]),]
  a_IYG_S=a[order(-a[,5]),]
  a_REG_S=a[order(-a[,6]),]
  a_bet=a[order(-a[,7]),]
  a_clos=a[order(-a[,8]),]
  a_klein=a[order(-a[,9]),]
  a_sim.jacard=a[order(-a[,10]),]
  a_sim.dic=a[order(-a[,11]),]
  a_sim.we=a[order(-a[,12]),]
  a_short=a[order(a[,13]),]
  a_pvalue_Top=a[order(a[,14]),]
  a_pvalue_first=a[order(a[,15]),]
  a_Jaccard_Q=a[order(-a[,16]),]
  a_Jaccard_N=a[order(-a[,17]),]
  a_Jaccard=a[order(-a[,18]),] 
  a_NB=a[order(-a[,20]),]
  a_RF=a[order(-a[,21]),]
  a_gbm=a[order(-a[,22]),]
  a_log=a[order(-a[,23]),]
  a_svm=a[order(-a[,24]),]
  a_lda=a[order(-a[,25]),]
  a_ens=a[order(-a[,26]),]
  
  
  # Rank of each left-out gene based on different criteria
  for (j in 1:length(cregiyg)){
	b_deg[j,i]=which(a_deg[,1]==cregiyg[j])
	b_IYG[j,i]=which(a_IYG[,1]==cregiyg[j])
	b_IYG_S[j,i]=which(a_IYG_S[,1]==cregiyg[j])
	b_REG[j,i]=which(a_REG[,1]==cregiyg[j])
	b_REG_S[j,i]=which(a_REG_S[,1]==cregiyg[j])
	b_bet[j,i]=which(a_bet[,1]==cregiyg[j])
	b_clos[j,i]=which(a_clos[,1]==cregiyg[j])
	b_klein[j,i]=which(a_klein[,1]==cregiyg[j])
	b_sim.jacard[j,i]=which(a_sim.jacard[,1]==cregiyg[j])
	b_sim.dic[j,i]=which(a_sim.dic[,1]==cregiyg[j])
	b_sim.we[j,i]=which(a_sim.we[,1]==cregiyg[j])
	b_short[j,i]=which(a_short[,1]==cregiyg[j])
	b_pvalue_Top[j,i]=which(a_pvalue_Top[,1]==cregiyg[j])
	b_pvalue_first[j,i]=which(a_pvalue_first[,1]==cregiyg[j])
	b_Jaccard_Q[j,i]=which(a_Jaccard_Q[,1]==cregiyg[j])
	b_Jaccard_N[j,i]=which(a_Jaccard_N[,1]==cregiyg[j])
	b_Jaccard[j,i]=which(a_Jaccard[,1]==cregiyg[j])
	b_NB[j,i]=which(a_NB[,1]==cregiyg[j])
	b_RF[j,i]=which(a_RF[,1]==cregiyg[j])
	b_gbm[j,i]=which(a_gbm[,1]==cregiyg[j])
	b_log[j,i]=which(a_log[,1]==cregiyg[j])
	b_svm[j,i]=which(a_svm[,1]==cregiyg[j])
	b_lda[j,i]=which(a_lda[,1]==cregiyg[j])
	b_ens[j,i]=which(a_ens[,1]==cregiyg[j])
  }
}


com=cbind(diag(b_deg),diag(b_IYG),diag(b_IYG_S),diag(b_REG),diag(b_REG_S),diag(b_bet),diag(b_clos),diag(b_klein),diag(b_sim.jacard)
,diag(b_sim.dic),diag(b_sim.we),diag(b_short),diag(b_pvalue_Top),diag(b_pvalue_first),diag(b_Jaccard_Q),diag(b_Jaccard_N),diag(b_Jaccard),
diag(b_NB),diag(b_RF),diag(b_gbm),diag(b_log),diag(b_svm),diag(b_lda),diag(b_ens))

colnames(com)=c("Degree","IYG","IYG_shared","REG","REG_shared","Betwee","clos","klein","sim.jacard","sim.dic","sim.we","short","pvalue_Top","pvalue_first",
"Jaccard_Q","Jaccard_N","Jaccard","NaiveBayes","RF","GBM","Logistic","SVM","LDA","Ensemble")

write.table(com,file="~/prioritization/Leave-one-out/Total/Results/results_Ens_GO.txt",sep="\t",col.names=T,row.names=F,quote=F)

