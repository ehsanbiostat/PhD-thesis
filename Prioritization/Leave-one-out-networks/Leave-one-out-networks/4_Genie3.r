load("~/prioritization/Combined_Network/New/Leave_one_out/Genie3.RData")
a1=a1_graph
cregiyg=read.table("~/prioritization/Combined_Network/New/Filtering_list/cregiyg.txt")
library(Matrix)
library(igraph)
sim.jacard=similarity.jaccard(a1)
sim.jacard=sim.jacard[,unlist(cregiyg[,2])]
save(sim.jacard,file="~/prioritization/Leave-one-out/New/Leave-one-out-network/4_Genie3.RData")

