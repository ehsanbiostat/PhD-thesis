a=read.table("~/prioritization/Leave-one-out/New/Total_topology/global.txt",header=T)
cregiyg=read.table("~/prioritization/Combined_Network/New/Filtering_list/cregiyg.txt") # Known genes list (IYG & REG)
b=read.table("~/prioritization/Leave-one-out/New/GO/GO_Total.txt",header=T)  # GO analysis
load("~/prioritization/Leave-one-out/New/GO/GO_Total_separate.RData")

# Upload packages
library(e1071)
library(randomForest)
library(gbm)
library(Matrix)
library(igraph)
library(MASS)
library(glmnet)

# Creating response vector
cregiyg=unlist(cregiyg[,2])
y=sparseMatrix(c(cregiyg,27376),rep(1,148),x=1)
y[27376,]=0
y=as.matrix(y)


# Creating ranking matrix
b_bet=matrix(NA,length(cregiyg),1)
b_deg=matrix(NA,length(cregiyg),1)
b_cregiyg=matrix(NA,length(cregiyg),1)
b_cregiyg_S=matrix(NA,length(cregiyg),1)
b_clos=matrix(NA,length(cregiyg),1)
b_klein=matrix(NA,length(cregiyg),1)
b_sim.jacard=matrix(NA,length(cregiyg),1)
b_sim.dic=matrix(NA,length(cregiyg),1)
b_sim.we=matrix(NA,length(cregiyg),1)
b_short=matrix(NA,length(cregiyg),1)
b_pvalue_Top=matrix(NA,length(cregiyg),1)
b_pvalue_first=matrix(NA,length(cregiyg),1)
b_Jaccard_Q=matrix(NA,length(cregiyg),1)
b_Jaccard_N=matrix(NA,length(cregiyg),1)
b_Jaccard=matrix(NA,length(cregiyg),1)
b_NB=matrix(NA,length(cregiyg),1)
b_RF=matrix(NA,length(cregiyg),1)
b_gbm=matrix(NA,length(cregiyg),1)
b_svm=matrix(NA,length(cregiyg),1)
b_lda=matrix(NA,length(cregiyg),1)
b_glmnet=matrix(NA,length(cregiyg),1)
b_ens=matrix(NA,length(cregiyg),1)
b_ens1=matrix(NA,length(cregiyg),1)
b_ens2=matrix(NA,length(cregiyg),1)
b_ens3=matrix(NA,length(cregiyg),1)
b_ens4=matrix(NA,length(cregiyg),1)
b_ens5=matrix(NA,length(cregiyg),1)
b_ens6=matrix(NA,length(cregiyg),1)
b_ens7=matrix(NA,length(cregiyg),1)
b_ens8=matrix(NA,length(cregiyg),1)
r=matrix(NA,27376,6)
b_ens1_R=matrix(NA,length(cregiyg),1)
b_ens2_R=matrix(NA,length(cregiyg),1)
b_ens3_R=matrix(NA,length(cregiyg),1)
b_ens4_R=matrix(NA,length(cregiyg),1)
b_ens5_R=matrix(NA,length(cregiyg),1)
b_ens6_R=matrix(NA,length(cregiyg),1)
b_ens7_R=matrix(NA,length(cregiyg),1)
b_ens8_R=matrix(NA,length(cregiyg),1)
b_ens9_R=matrix(NA,length(cregiyg),1)
b_ens10_R=matrix(NA,length(cregiyg),1)


# Reading each left-out dataset
  s=data.frame(w[1],w[2],w[3])
  s=s[,-c(2,4,6,11,13,15,20,22,24)]
  colnames(s)=c("PTC","PFC","P5C","JQC","JNC","JC","PTF","PFF","P5F","JQF","JNF","JF","PTP","PFP","P5P","JQP","JNP","JP")
  a=cbind(c(1:dim(a)[1]),a,b[,-c(2,4,6)],s,y) # Combination of data sets except odd ratio values
  a[a[,11]=="Inf",11]=1000
  a1=a
    
  # Datasets
  #new=test=data.frame(a[cregiyg[i],-c(1,dim(a)[2])])
  x=data.frame(a1[,-c(1,dim(a1)[2])])
  Y=data.frame(a1[,dim(a)[2]])
  
  # Naive Bayes
  NB=naiveBayes(x,unlist(Y))
  #p=predict(NB,new,type=c("raw"))[2]
  fit=predict(NB,x,type=c("raw"))[,2]
  a=cbind(a,fit)
  colnames(a)[dim(a)[2]]=c("NaiveBayes")
    
  # Random Forest
  rf=randomForest(x,as.factor(unlist(Y)),ntree=1000)
  fit=rf$votes[,2]
  #p=predict(rf,new,type="vote")[2]
  a=cbind(a,fit)
  colnames(a)[dim(a)[2]]=c("RF")
  
  # GBM
  gb=gbm.fit(x,unlist(Y),distribution="adaboost",n.trees=100)
  #p=predict(gb,new,n.trees=100)
  fit=gb$fit
  a=cbind(a,fit)
  colnames(a)[dim(a)[2]]=c("GBM")
  
     
  # SVM
  sv=svm(x,Y,probability=TRUE)
  fit=sv$fitted
  #p=predict(sv,new,type="response")
  a=cbind(a,fit)
  colnames(a)[dim(a)[2]]=c("SVM")
  
  # Linear Discriminante Analysis (Fisher)
  LDA=lda(x[,-5],unlist(Y)) # Since there is no adequate variation in clos feature LDA cannot use it
  fit=predict(LDA,method=c("predictive"))$posterior[,2]
  #p=predict(LDA,new[,-5],method=c("predictive"))$posterior[,2]
  a=cbind(a,fit)
  colnames(a)[dim(a)[2]]=c("LDA")
  
  
  # GLMNET
  glmn=cv.glmnet(as.matrix(x),as.factor(as.vector(unlist(Y))),family=c("binomial"),nlambda=90)
  fit=predict(glmn, as.matrix(x),type=c("response"))
  #p=predict(glmn, as.matrix(new),type=c("response"))
  a=cbind(a,fit)
  colnames(a)[dim(a)[2]]=c("Glmnet")
  
  
  # Ensemble
  Ens=a$RF+a$SVM+a$LDA
  a=cbind(a,Ens)
  
  
  # Ensemble 1
  Ens1=a$RF+a$LDA
  a=cbind(a,Ens1)
  
  # Ensemble 2
  Ens2=a$RF+a$NaiveBayes+a$LDA
  a=cbind(a,Ens2)
  
  # Ensemble 3
  Ens3=a$RF+a$SVM+a$Glmnet
  a=cbind(a,Ens3)
  
  # Ensemble 4
  Ens4=a$RF+a$SVM+a$LDA+a$NaiveBayes
  a=cbind(a,Ens4)
  
  # Ensemble 5
  Ens5=a$RF+a$SVM+a$LDA+a$NaiveBayes+a$Glmnet
  a=cbind(a,Ens5)
  
  # Ensemble 6
  Ens6=a$RF+a$LDA+a$Glmnet
  a=cbind(a,Ens6)
  
  # Ensemble 7
  Ens7=a$RF+a$Glmnet
  a=cbind(a,Ens7)
  
  # Ensemble 8
  Ens8=a$LDA+a$Glmnet
  a=cbind(a,Ens8)
  
  # Ranking based on different criteria
  a_deg=a[order(-a[,2]),]
  a_cregiyg=a[order(-a[,3]),]
  a_cregiyg_S=a[order(-a[,4]),]
  a_bet=a[order(-a[,5]),]
  a_clos=a[order(-a[,6]),]
  a_klein=a[order(-a[,7]),]
  a_sim.jacard=a[order(-a[,8]),]
  a_sim.dic=a[order(-a[,9]),]
  a_sim.we=a[order(-a[,10]),]
  a_short=a[order(a[,11]),]
  a_pvalue_Top=a[order(a[,12]),]
  a_pvalue_first=a[order(a[,13]),]
  a_pvalue_5=a[order(a[,14]),]
  a_Jaccard_Q=a[order(-a[,15]),]
  a_Jaccard_N=a[order(-a[,16]),]
  a_Jaccard=a[order(-a[,17]),] 
  a_NB=a[order(-a[,37]),]
  a_RF=a[order(-a[,38]),]
  a_gbm=a[order(-a[,39]),]
  a_svm=a[order(-a[,40]),]
  a_lda=a[order(-a[,41]),]
  a_glmnet=a[order(-a[,42]),]
  a_ens=a[order(-a[,43]),]
  a_ens1=a[order(-a[,44]),]
  a_ens2=a[order(-a[,45]),]
  a_ens3=a[order(-a[,46]),]
  a_ens4=a[order(-a[,47]),]
  a_ens5=a[order(-a[,48]),]
  a_ens6=a[order(-a[,49]),]
  a_ens7=a[order(-a[,50]),]
  a_ens8=a[order(-a[,51]),]
  ens=cbind(a_NB[,1],a_RF[,1],a_gbm[,1],a_svm[,1],a_lda[,1],a_glmnet[,1])
  
  for (k in 1:6){
    for (j in 1:27376){
	 r[j,k]=which(j==ens[,k])
	}
  }
 colnames(r)=c("NB","RF","GBM","SVM","LDA","Glmnet")
 r=as.data.frame(r)
 ens1_R=r$NB+r$RF+r$GBM+r$SVM+r$LDA+r$Glmnet
 ens2_R=r$NB+r$RF
 ens3_R=r$RF+r$LDA+r$Glmnet
 ens4_R=r$RF+r$LDA
 ens5_R=r$RF+r$Glmnet
 ens6_R=r$LDA+r$Glmnet
 ens7_R=r$RF+r$SVM+r$LDA+r$Glmnet
 ens8_R=r$NB+r$RF+r$SVM
 ens9_R=r$NB+r$RF+r$GBM+r$SVM
 ens10_R=r$NB+r$RF+r$LDA+r$Glmnet
 ens_total=cbind(ens1_R,ens2_R,ens3_R,ens4_R,ens5_R,ens6_R,ens7_R,ens8_R,ens9_R,ens10_R,a[,1])
 a=cbind(a,ens_total[,-11])
 
   i=1
  # Rank of each left-out gene based on different criteria
  for (j in 1:length(cregiyg)){
	b_deg[j,i]=which(a_deg[,1]==cregiyg[j])
	b_cregiyg[j,i]=which(a_cregiyg[,1]==cregiyg[j])
	b_cregiyg_S[j,i]=which(a_cregiyg_S[,1]==cregiyg[j])
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
	b_svm[j,i]=which(a_svm[,1]==cregiyg[j])
	b_lda[j,i]=which(a_lda[,1]==cregiyg[j])
	b_glmnet[j,i]=which(a_glmnet[,1]==cregiyg[j])
	b_ens[j,i]=which(a_ens[,1]==cregiyg[j])
	b_ens1[j,i]=which(a_ens1[,1]==cregiyg[j])
	b_ens2[j,i]=which(a_ens2[,1]==cregiyg[j])
	b_ens3[j,i]=which(a_ens3[,1]==cregiyg[j])
	b_ens4[j,i]=which(a_ens4[,1]==cregiyg[j])
	b_ens5[j,i]=which(a_ens5[,1]==cregiyg[j])
	b_ens6[j,i]=which(a_ens6[,1]==cregiyg[j])
	b_ens7[j,i]=which(a_ens7[,1]==cregiyg[j])
	b_ens8[j,i]=which(a_ens8[,1]==cregiyg[j])
	b_ens1_R[j,i]=which(ens_total[order(ens_total[,1]),dim(ens_total)[2]]==cregiyg[j])
	b_ens2_R[j,i]=which(ens_total[order(ens_total[,2]),dim(ens_total)[2]]==cregiyg[j])
	b_ens3_R[j,i]=which(ens_total[order(ens_total[,3]),dim(ens_total)[2]]==cregiyg[j])
	b_ens4_R[j,i]=which(ens_total[order(ens_total[,4]),dim(ens_total)[2]]==cregiyg[j])
	b_ens5_R[j,i]=which(ens_total[order(ens_total[,5]),dim(ens_total)[2]]==cregiyg[j])
	b_ens6_R[j,i]=which(ens_total[order(ens_total[,6]),dim(ens_total)[2]]==cregiyg[j])
	b_ens7_R[j,i]=which(ens_total[order(ens_total[,7]),dim(ens_total)[2]]==cregiyg[j])
	b_ens8_R[j,i]=which(ens_total[order(ens_total[,8]),dim(ens_total)[2]]==cregiyg[j])
	b_ens9_R[j,i]=which(ens_total[order(ens_total[,9]),dim(ens_total)[2]]==cregiyg[j])
	b_ens10_R[j,i]=which(ens_total[order(ens_total[,10]),dim(ens_total)[2]]==cregiyg[j])
  }

  
com=data.frame(b_bet,b_deg,b_cregiyg,b_cregiyg_S,b_clos,b_klein,b_sim.jacard,b_sim.dic,b_sim.we,b_short,b_pvalue_Top,b_pvalue_first,b_Jaccard_Q,b_Jaccard_N,b_Jaccard,b_NB,b_RF,b_gbm,b_svm,b_lda,b_glmnet,b_ens,b_ens1,b_ens2,b_ens3,b_ens4,b_ens5,b_ens6,b_ens7,b_ens8,b_ens1_R,b_ens2_R,b_ens3_R,b_ens4_R,b_ens5_R,b_ens6_R,b_ens7_R,b_ens8_R,b_ens9_R,b_ens10_R)

colnames(com)=c("Degree","cregiyg","vregiyg_shared","Betwee","clos","klein","sim.jacard","sim.dic","sim.we","short","pvalue_Top","pvalue_first",
"Jaccard_Q","Jaccard_N","Jaccard","NaiveBayes","RF","GBM","SVM","LDA","Glmnet","RF+SVM+LDA","RF+LDA","RF+NaiveBayes+LDA","RF+SVM+Glmnet","RF+SVM+LDA+NaiveBayes","RF+SVM+LDA+NaiveBayes+Glmnet","RF+LDA+Glmnet","RF+Glmnet","LDA+Glmnet",
"ens1","ens2","ens3","ens4","ens5","ens6","ens7","ens8","ens9","ens10")


write.table(com,file="~/prioritization/Ranking/GR_Ranking.txt",sep="\t",col.names=T,row.names=F)
write.table(a,file="~/prioritization/Ranking/Total.txt",sep="\t",col.names=T,row.names=F)