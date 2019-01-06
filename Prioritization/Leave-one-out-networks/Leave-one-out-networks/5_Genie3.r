load("~/prioritization/Combined_Network/New/Leave_one_out/Genie3.RData")
a1=a1_graph
cregiyg=read.table("~/prioritization/Combined_Network/New/Filtering_list/cregiyg.txt")
library(Matrix)
library(igraph)
sim.dic=similarity.dice(a1)
sim.dic=sim.dic[,unlist(cregiyg[,2])]
save(sim.dic,file="~/prioritization/Leave-one-out/New/Leave-one-out-network/5_Genie3.RData")

