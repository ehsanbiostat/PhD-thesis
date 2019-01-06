#!/usr/local/bin/perl -w

use strict;

for(my $i=1;$i<148;$i++){
  my $jobfile="PPI_Leave".$i.".r";
  open(OUT,">$jobfile");
  print OUT  "load(\"~/prioritization/Combined_Network/New/Leave_one_out/PPI.RData\")
cregiyg=read.table(\"~/prioritization/Combined_Network/New/Filtering_list/cregiyg.txt\")
library(Matrix)
library(igraph)
a=get.adjacency(a1_graph,type=c(\"both\"),sparse=T)
a1=a1_graph
cregiyg=unlist(cregiyg[,2])[-".$i."]
p=a[,unlist(cregiyg)]
d=p%*%t(p)    #Matrix multiplication shared genes with cregiyg
sh=d[,unlist(cregiyg)]
cregiyg_sh=apply(sh,1,sum)
cregiyg=a[,unlist(cregiyg)]
decregiyg=apply(cregiyg,1,sum) #Number of connection with cregiyg

f=data.frame(decregiyg,cregiyg_sh)
write.table(f,file=\"~/prioritization/Leave-one-out/New/Leave-one-out-network/0_PPI".$i.".txt\",sep=\"\\t\",col.names=T,row.names=F)";
   close OUT;
}

