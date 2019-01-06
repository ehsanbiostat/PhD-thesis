load("~/prioritization/Combined_Network/New/Total/Total_matrix.RData")
load("~/prioritization/Combined_Network/New/Total/Total_graph.RData")
cregiyg=read.table("~/prioritization/Combined_Network/New/Filtering_list/cregiyg.txt")
cregiyg=cregiyg[,2]
a=Total_matrix
a1=Total_graph
library(Matrix)
library(igraph)
cregiyg=as.matrix(cregiyg)
bet=betweenness(a1,directed=F)  # Betweennees Centrality

save(bet,file="~/prioritization/Leave-one-out/New/Total_topology/Bet.RData")