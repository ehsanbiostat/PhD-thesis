load("~/prioritization/Combined_Network/New/Leave_one_out/AG.RData")
a1=a1_graph
cregiyg=read.table("~/prioritization/Combined_Network/New/Filtering_list/cregiyg.txt")
library(Matrix)
library(igraph)
short=shortest.paths(a1,c(1:length(a1_graph[1])),c(unlist(cregiyg[,2])))
save(short,file="~/prioritization/Leave-one-out/New/Leave-one-out-network/7_AG.RData")

