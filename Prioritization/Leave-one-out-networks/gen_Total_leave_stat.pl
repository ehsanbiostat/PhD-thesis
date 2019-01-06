#!/usr/local/bin/perl -w

use strict;

for(my $i=57;$i<147;$i++){
  my $jobfile="Total_Leave".$i.".r";
  open(OUT,">$jobfile");
  print OUT  "load(\"~/../../group/biocomp/users/ehsab/Combined_Network/Leave-one-out/Total_Genemania.RData\")
iyg=read.table(\"~/prioritization/Combined_Network/Filtering_list/IYG.txt\")
reg=read.table(\"~/prioritization/Combined_Network/Filtering_list/REG.txt\")
cregiyg=read.table(\"~/prioritization/Combined_Network/Filtering_list/cregiyg.txt\")
a=Total_AG_matrix
a1=Total_AG_graph
#iyg=unlist(iyg)[-".$i."]
reg=unlist(reg)[-".($i-56)."]
cregiyg=unlist(cregiyg)[-".$i."]
library(Matrix)
library(igraph)
p=a[,unlist(iyg)]
d=p%*%t(p)   #Matrix multiplication shared genes with IYG
sh=d[,unlist(iyg)]
IYG_S=apply(sh,1,sum)
p=a[,unlist(reg)]
d=p%*%t(p)    #Matrix multiplication shared genes with REG
sh=d[,unlist(reg)]
REG_S=apply(sh,1,sum)
p=a[,unlist(cregiyg)]
d=p%*%t(p)    #Matrix multiplication shared genes with cregiyg
sh=d[,unlist(cregiyg)]
CREGIYG_S=apply(sh,1,sum)
ciyg=a[,unlist(iyg)]
creg=a[,unlist(reg)]
cregiyg=a[,unlist(cregiyg)]
deiyg=apply(ciyg,1,sum) #Number of connection with IYG
dereg=apply(creg,1,sum) #Number of connection with REG
decregiyg=apply(cregiyg,1,sum) #Number of connection with cregiyg

f=cbind(as.matrix(deiyg),as.matrix(dereg),as.matrix(IYG_S),as.matrix(REG_S),as.matrix(decregiyg),as.matrix(CREGIYG_S))

colnames(f)=c(\"IYG\",\"REG\",\"shared_IYG\",\"shared_REG\",\"cregiyg\",\"shared_cregiyg\")
write.table(f,file=\"~/prioritization/Leave-one-out/Total/Leave-one-out-network/Genemania/re".$i.".txt\",sep=\"\\t\",col.names=T,row.names=F)";
   close OUT;
}

