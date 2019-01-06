load("~/prioritization/Combined_Network/New/Total/Total_matrix.RData")
load("~/prioritization/Combined_Network/New/Total/Total_graph.RData")
cregiyg=read.table("~/prioritization/Combined_Network/New/Filtering_list/cregiyg.txt")
cregiyg=cregiyg[,2]
a=Total_matrix
a1=Total_graph
library(Matrix)
library(igraph)
clos=closeness(a1,mode=c("all")) #Closeness 

save(clos,file="~/prioritization/Leave-one-out/New/Total_topology/closeness.RData")