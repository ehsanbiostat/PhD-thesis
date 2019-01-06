#!/usr/local/bin/perl -w

use strict;

for(my $i=1;$i<148;$i++){
  my $jobfile="GO".$i."_text_Leave.r";
  open(OUT,">$jobfile");
  print OUT  "load(\"~/prioritization/Leave-one-out/New/Leave-one-out-network/text.RData\")  # Whole Network data sets
d=w
cregiyg=read.table(\"~/prioritization/Combined_Network/New/Filtering_list/cregiyg.txt\") # Known genes list (IYG & REG)
load(\"~/../../group/biocomp/users/ehsab/Leave-one-out/GO/text_Leave-one-out.RData\")  # GO analysis
load(\"~/../../group/biocomp/users/ehsab/Leave-one-out/GO/text_GO_separated.RData\")

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
b_bet=matrix(NA,length(cregiyg),length(d))
b_deg=matrix(NA,length(cregiyg),length(d))
b_cregiyg=matrix(NA,length(cregiyg),length(d))
b_cregiyg_S=matrix(NA,length(cregiyg),length(d))
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
b_log1=matrix(NA,length(cregiyg),length(d))
b_svm=matrix(NA,length(cregiyg),length(d))
b_lda=matrix(NA,length(cregiyg),length(d))
b_glmnet=matrix(NA,length(cregiyg),length(d))
b_ens=matrix(NA,length(cregiyg),length(d))
b_ens1=matrix(NA,length(cregiyg),length(d))
b_ens2=matrix(NA,length(cregiyg),length(d))
b_ens3=matrix(NA,length(cregiyg),length(d))
b_ens4=matrix(NA,length(cregiyg),length(d))
b_ens5=matrix(NA,length(cregiyg),length(d))
b_ens6=matrix(NA,length(cregiyg),length(d))
b_ens7=matrix(NA,length(cregiyg),length(d))
b_ens8=matrix(NA,length(cregiyg),length(d))
r=matrix(NA,27376,8)
b_ens1_R=matrix(NA,length(cregiyg),length(d))
b_ens2_R=matrix(NA,length(cregiyg),length(d))
b_ens3_R=matrix(NA,length(cregiyg),length(d))
b_ens4_R=matrix(NA,length(cregiyg),length(d))
b_ens5_R=matrix(NA,length(cregiyg),length(d))
b_ens6_R=matrix(NA,length(cregiyg),length(d))
b_ens7_R=matrix(NA,length(cregiyg),length(d))
b_ens8_R=matrix(NA,length(cregiyg),length(d))
b_ens9_R=matrix(NA,length(cregiyg),length(d))
b_ens10_R=matrix(NA,length(cregiyg),length(d))


# Reading each left-out dataset
  i=".$i."
  print(i)
  a=as.data.frame(d[i])
  b=as.data.frame(g[i])
  s=as.data.frame(sep[i])
  s=s[,-c(2,4,6,11,13,15,20,22,24)]
  colnames(s)=c(\"PTC\",\"PFC\",\"P5C\",\"JQC\",\"JNC\",\"JC\",\"PTF\",\"PFF\",\"P5F\",\"JQF\",\"JNF\",\"JF\",\"PTP\",\"PFP\",\"P5P\",\"JQP\",\"JNP\",\"JP\")
  a=cbind(c(1:dim(a)[1]),a,b[,-c(2,4,6)],s,y) # Combination of two data sets except odd ratio values
  a[a[,11]==\"Inf\",11]=1000
  a1=a[-cregiyg[i],]
    
  # Datasets
  new=test=data.frame(a[cregiyg[i],-c(1,dim(a)[2])])
  x=data.frame(a1[,-c(1,dim(a1)[2])])
  Y=data.frame(a1[,dim(a)[2]])
  
  # Naive Bayes
  NB=naiveBayes(x,unlist(Y))
  p=predict(NB,new,type=c(\"raw\"))[2]
  fit=predict(NB,x,type=c(\"raw\"))[,2]
  a=cbind(a,c(fit[1:c(cregiyg[i]-1)],p,fit[cregiyg[i]:length(fit)]))
  colnames(a)[dim(a)[2]]=c(\"NaiveBayes\")
    
  # Random Forest
  rf=randomForest(x,as.factor(unlist(Y)),ntree=1000)
  fit=rf\$votes[,2]
  p=predict(rf,new,type=\"vote\")[2]
  a=cbind(a,c(fit[1:c(cregiyg[i]-1)],p,fit[cregiyg[i]:length(fit)]))
  colnames(a)[dim(a)[2]]=c(\"RF\")
  
  # GBM
  gb=gbm.fit(x,unlist(Y),distribution=\"adaboost\",n.trees=100)
  p=predict(gb,new,n.trees=100)
  fit=gb\$fit
  a=cbind(a,c(fit[1:c(cregiyg[i]-1)],p,fit[cregiyg[i]:length(fit)]))
  colnames(a)[dim(a)[2]]=c(\"GBM\")
  
     
  # Logistic Regression 
  data=data.frame(Y=a1\$y,alldeg=a1[,2],decregiyg=a1[,3],cregiyg_sh=a1[,4],bet=a1[,5],clos=a1[,6]
  ,klein=a1[,7],sim.jacard=a1[,8],sim.dic=a1[,9],sim.we=a1[,10],short=a1[,11],pvalue_TopTen=a1[,12],pvalue_first=a1[,13],pvalue_5=a1[,14],Jaccard_Q=a1[,15],Jaccard_N=a1[,16],Jaccard=a1[,17]
  ,PTC=a1[,18],PFC=a1[,19],P5C=a1[,20],JQC=a1[,21],JNC=a1[,22],JC=a1[,23],PTF=a1[,24],PFF=a1[,25],P5F=a1[,26],JQF=a1[,27],JNF=a1[,28]
  ,JF=a1[,29],PTP=a1[,30],PFP=a1[,35],P5P=a1[,32],JQP=a1[,33],JNP=a1[,34],JP=a1[,35])
  model=glm(Y~alldeg+decregiyg+cregiyg_sh+bet+clos+klein+sim.jacard+sim.dic+sim.we+short
  +pvalue_TopTen+pvalue_first+pvalue_5+Jaccard_Q+Jaccard_N+Jaccard+PTC+PFC+P5C+JQC+JNC+JC+PTF+PFF+P5F+JQF+JNF+JF+PTP+PFP+P5P+JQP+JNP+JP,data=data,family=\"binomial\")
  newglm=a[cregiyg[i],c(2:35)]
  p=predict(model,newglm,type=c(\"response\"))
  fit=predict(model,type=c(\"response\"))
  a=cbind(a,c(fit[1:c(cregiyg[i]-1)],p,fit[cregiyg[i]:length(fit)]))
  colnames(a)[dim(a)[2]]=c(\"Logistic\")
  
  
  # Logistic Regression 1
  data=data.frame(Y=a1\$y,alldeg=a1[,2],decregiyg=a1[,3],cregiyg_sh=a1[,4],bet=a1[,5],clos=a1[,6]
  ,klein=a1[,7],sim.jacard=a1[,8],sim.dic=a1[,9],sim.we=a1[,10],short=a1[,11],PTC=a1[,18],PFC=a1[,19],P5C=a1[,20],JQC=a1[,21],JNC=a1[,22],JC=a1[,23],PTF=a1[,24],PFF=a1[,25],P5F=a1[,26],JQF=a1[,27],JNF=a1[,28]
  ,JF=a1[,29],PTP=a1[,30],PFP=a1[,35],P5P=a1[,32],JQP=a1[,33],JNP=a1[,34],JP=a1[,35])
  model=glm(Y~alldeg+decregiyg+cregiyg_sh+bet+clos+klein+sim.jacard+sim.dic+sim.we+short
  +PTC+PFC+P5C+JQC+JNC+JC+PTF+PFF+P5F+JQF+JNF+JF+PTP+PFP+P5P+JQP+JNP+JP,data=data,family=\"binomial\")
  newglm=a[cregiyg[i],c(2:11,18:35)]
  p=predict(model,newglm,type=c(\"response\"))
  fit=predict(model,type=c(\"response\"))
  a=cbind(a,c(fit[1:c(cregiyg[i]-1)],p,fit[cregiyg[i]:length(fit)]))
  colnames(a)[dim(a)[2]]=c(\"Logistic1\")
  
  
  
  # SVM
  sv=svm(x,Y,probability=TRUE)
  fit=sv\$fitted
  p=predict(sv,new,type=\"response\")
  a=cbind(a,c(fit[1:c(cregiyg[i]-1)],p,fit[cregiyg[i]:length(fit)]))
  colnames(a)[dim(a)[2]]=c(\"SVM\")
  
  # Linear Discriminante Analysis (Fisher)
  LDA=lda(x[,-5],unlist(Y))
  fit=predict(LDA,method=c(\"predictive\"))\$posterior[,2]
  p=predict(LDA,new[,-5],method=c(\"predictive\"))\$posterior[,2]
  a=cbind(a,c(fit[1:c(cregiyg[i]-1)],p,fit[cregiyg[i]:length(fit)]))
  colnames(a)[dim(a)[2]]=c(\"LDA\")
  
  
  # GLMNET
  glmn=cv.glmnet(as.matrix(x),as.factor(as.vector(unlist(Y))),family=c(\"binomial\"),nlambda=90)
  fit=predict(glmn, as.matrix(x),type=c(\"response\"))
  p=predict(glmn, as.matrix(new),type=c(\"response\"))
  a=cbind(a,c(fit[1:c(cregiyg[i]-1)],p,fit[cregiyg[i]:length(fit)]))
  colnames(a)[dim(a)[2]]=c(\"Glmnet\")
  
  
  # Ensemble
  Ens=a\$RF+a\$Logistic+a\$LDA
  a=cbind(a,Ens)
  
  
  # Ensemble 1
  Ens1=a\$RF+a\$LDA
  a=cbind(a,Ens1)
  
  # Ensemble 2
  Ens2=a\$RF+a\$NaiveBayes+a\$LDA
  a=cbind(a,Ens2)
  
  # Ensemble 3
  Ens3=a\$RF+a\$Logistic1+a\$LDA
  a=cbind(a,Ens3)
  
  # Ensemble 4
  Ens4=a\$RF+a\$Logistic1+a\$LDA+a\$NaiveBayes
  a=cbind(a,Ens4)
  
  # Ensemble 5
  Ens5=a\$RF+a\$Logistic1+a\$LDA+a\$NaiveBayes+a\$Glmnet
  a=cbind(a,Ens5)
  
  # Ensemble 6
  Ens6=a\$RF+a\$LDA+a\$Glmnet
  a=cbind(a,Ens6)
  
  # Ensemble 7
  Ens7=a\$RF+a\$Glmnet
  a=cbind(a,Ens7)
  
  # Ensemble 8
  Ens8=a\$LDA+a\$Glmnet
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
  a_log=a[order(-a[,40]),]
  a_log1=a[order(-a[,41]),]
  a_svm=a[order(-a[,42]),]
  a_lda=a[order(-a[,43]),]
  a_glmnet=a[order(-a[,44]),]
  a_ens=a[order(-a[,45]),]
  a_ens1=a[order(-a[,46]),]
  a_ens2=a[order(-a[,47]),]
  a_ens3=a[order(-a[,48]),]
  a_ens4=a[order(-a[,49]),]
  a_ens5=a[order(-a[,50]),]
  a_ens6=a[order(-a[,51]),]
  a_ens7=a[order(-a[,52]),]
  a_ens8=a[order(-a[,53]),]
  ens=cbind(a_NB[,1],a_RF[,1],a_gbm[,1],a_log[,1],a_log1[,1],a_svm[,1],a_lda[,1],a_glmnet[,1])
  
  for (k in 1:8){
    for (j in 1:27376){
	 r[j,k]=which(j==ens[,k])
	}
  }
 colnames(r)=c(\"NB\",\"RF\",\"GBM\",\"log\",\"log_Q\",\"SVM\",\"LDA\",\"Glmnet\")
 r=as.data.frame(r)
 ens1_R=r\$NB+r\$RF+r\$GBM+r\$log+r\$log_Q+r\$SVM+r\$LDA+r\$Glmnet
 ens2_R=r\$NB+r\$RF
 ens3_R=r\$RF+r\$LDA+r\$Glmnet
 ens4_R=r\$RF+r\$LDA
 ens5_R=r\$RF+r\$Glmnet
 ens6_R=r\$LDA+r\$Glmnet
 ens7_R=r\$RF+r\$log_Q+r\$LDA+r\$Glmnet
 ens8_R=r\$NB+r\$RF+r\$log
 ens9_R=r\$NB+r\$RF+r\$GBM+r\$log
 ens10_R=r\$NB+r\$RF+r\$LDA+r\$Glmnet
 ens_total=cbind(ens1_R,ens2_R,ens3_R,ens4_R,ens5_R,ens6_R,ens7_R,ens8_R,ens9_R,ens10_R,a[,1])
 
   
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
	b_log[j,i]=which(a_log[,1]==cregiyg[j])
	b_log1[j,i]=which(a_log1[,1]==cregiyg[j])
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




com=cbind(diag(b_deg),diag(b_cregiyg),diag(b_cregiyg_S),diag(b_bet),diag(b_clos),diag(b_klein),diag(b_sim.jacard)
,diag(b_sim.dic),diag(b_sim.we),diag(b_short),diag(b_pvalue_Top),diag(b_pvalue_first),diag(b_Jaccard_Q),diag(b_Jaccard_N),diag(b_Jaccard),
diag(b_NB),diag(b_RF),diag(b_gbm),diag(b_log),diag(b_log1),diag(b_svm),diag(b_lda),diag(b_glmnet),diag(b_ens),diag(b_ens1),diag(b_ens2),diag(b_ens3),diag(b_ens4),diag(b_ens5),diag(b_ens6),diag(b_ens7),diag(b_ens8),
diag(b_ens1_R),diag(b_ens2_R),diag(b_ens3_R),diag(b_ens4_R),diag(b_ens5_R),diag(b_ens6_R),diag(b_ens7_R),diag(b_ens8_R),diag(b_ens9_R),diag(b_ens10_R))

colnames(com)=c(\"Degree\",\"cregiyg\",\"vregiyg_shared\",\"Betwee\",\"clos\",\"klein\",\"sim.jacard\",\"sim.dic\",\"sim.we\",\"short\",\"pvalue_Top\",\"pvalue_first\",
\"Jaccard_Q\",\"Jaccard_N\",\"Jaccard\",\"NaiveBayes\",\"RF\",\"GBM\",\"Logistic\",\"Logistic_Jac\",\"SVM\",\"LDA\",\"Glmnet\",\"RF+Logistic+LDA\",\"RF+LDA\",\"RF+NaiveBayes+LDA\",\"RF+Logistic1+LDA\",\"RF+Logistic1+LDA+NaiveBayes\",\"RF+Logistic1+LDA+NaiveBayes+Glmnet\",\"RF+LDA+Glmnet\",\"RF+Glmnet\",\"LDA+Glmnet\",
\"ens1\",\"ens2\",\"ens3\",\"ens4\",\"ens5\",\"ens6\",\"ens7\",\"ens8\",\"ens9\",\"ens10\")


write.table(com,file=\"~/prioritization/Leave-one-out/New/Leave-one-out-network/Results/Genie3/text_".$i.".txt\",sep=\"\\t\",col.names=T,row.names=F)";
   close OUT;
}

