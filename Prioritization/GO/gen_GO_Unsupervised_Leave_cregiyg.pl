#!/usr/local/bin/perl -w

use strict;

for(my $i=1;$i<148;$i++){
  my $jobfile="GO".$i."_Unsupervised_Leave.r";
  open(OUT,">$jobfile");
  print OUT  "load(\"~/prioritization/Leave-one-out/New/Total_topology/Leave-one-out.RData\")  # Whole Network data sets
cregiyg=read.table(\"~/prioritization/Combined_Network/New/Filtering_list/cregiyg.txt\") # Known genes list (IYG & REG)

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
r=matrix(NA,27376,10)
b_ens1_R=matrix(NA,length(cregiyg),length(d))


# Reading each left-out dataset
  i=".$i."
  print(i)
  a=as.data.frame(d[i])
  a=cbind(c(1:dim(a)[1]),a,y) # Combination of two data sets except odd ratio values
  a[a[,11]==\"Inf\",11]=1000
  
    
  
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
  ens=cbind(a_deg[,1],a_cregiyg[,1],a_cregiyg_S[,1],a_bet[,1],a_clos[,1],a_klein[,1],a_sim.jacard[,1],a_sim.dic[,1],a_sim.we[,1],a_short[,1])
  
  for (k in 1:10){
    for (j in 1:27376){
	 r[j,k]=which(j==ens[,k])
	}
  }
 colnames(r)=c(\"deg\",\"cregiyg\",\"cregiyg_S\",\"bet\",\"clos\",\"klein\",\"sim.jacard\",\"sim.dic\",\"sim.we\",\"shortest\")
 r=as.data.frame(r)
 ens_total=cbind(r\$deg+r\$cregiyg+r\$cregiyg_S+r\$bet+r\$clos+r\$klein+r\$sim.jacard+r\$sim.dic+r\$sim.we+r\$shortest,a[,1])
 
   
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
	b_ens1_R[j,i]=which(ens_total[order(ens_total[,1]),dim(ens_total)[2]]==cregiyg[j])
  }




com=cbind(diag(b_deg),diag(b_cregiyg),diag(b_cregiyg_S),diag(b_bet),diag(b_clos),diag(b_klein),diag(b_sim.jacard)
,diag(b_sim.dic),diag(b_sim.we),diag(b_short),diag(b_ens1_R))

colnames(com)=c(\"Degree\",\"cregiyg\",\"vregiyg_shared\",\"Betwee\",\"clos\",\"klein\",\"sim.jacard\",\"sim.dic\",\"sim.we\",\"short\",\"ens1\")


write.table(com,file=\"~/prioritization/Leave-one-out/New/GO/Results/Unsu_ensemble".$i.".txt\",sep=\"\\t\",col.names=T,row.names=F)";
   close OUT;
}

