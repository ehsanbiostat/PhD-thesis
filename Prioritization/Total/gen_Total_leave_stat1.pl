#!/usr/local/bin/perl -w

use strict;

for(my $i=1;$i<147;$i++){
  my $jobfile="Total_Leave".$i.".r";
  open(OUT,">$jobfile");
  print OUT  "load(\"~/prioritization/Combined_Network/Total/TOTAL1.RData\")
cregiyg=read.table(\"~/prioritization/Combined_Network/Filtering_list/cregiyg.txt\")
a=Total_matrix
a1=Total_graph
cregiyg=unlist(cregiyg)[-".$i."]
library(Matrix)
library(igraph)
p=a[,unlist(cregiyg)]
d=p%*%t(p)    #Matrix multiplication shared genes with cregiyg
sh=d[,unlist(cregiyg)]
CREGIYG_S=apply(sh,1,sum)
cregiyg=a[,unlist(cregiyg)]
decregiyg=apply(cregiyg,1,sum) #Number of connection with cregiyg

f=cbind(as.matrix(decregiyg),as.matrix(CREGIYG_S))

colnames(f)=c(\"cregiyg\",\"shared_cregiyg\")
write.table(f,file=\"~/prioritization/Leave-one-out/Total/re".$i.".txt\",sep=\"\\t\",col.names=T,row.names=F)";
   close OUT;
}

