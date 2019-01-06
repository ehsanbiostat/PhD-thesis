load("~/prioritization/Combined_Network/New/Leave_one_out/Mutant.RData")
a1=a1_graph
cregiyg=read.table("~/prioritization/Combined_Network/New/Filtering_list/cregiyg.txt")
library(Matrix)
library(igraph)
sim.we=similarity.invlogweighted(a1)
sim.we=sim.we[,unlist(cregiyg[,2])]
save(sim.we,file="~/prioritization/Leave-one-out/New/Leave-one-out-network/6_Mutant.RData")

