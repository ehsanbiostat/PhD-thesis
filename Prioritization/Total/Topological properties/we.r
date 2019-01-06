load("/ngsprojects/iwt/ehsab/Prioritization/prioritization/Combined_Network/New/Total/Total_matrix.RData")
load("/ngsprojects/iwt/ehsab/Prioritization/prioritization/Combined_Network/New/Total/Total_graph.RData")
cregiyg = read.table("/ngsprojects/iwt/ehsab/Prioritization/Plant_Cell/Cross_validation/10_fold/Seed_genes/10.txt")
cregiyg=cregiyg[, 2]
a = Total_matrix
a1 = Total_graph
library(Matrix)
library(igraph)

sim.we=similarity.invlogweighted(a1)
sim.we=sim.we[,unlist(cregiyg)]

save(sim.we,file="/ngsprojects/iwt/ehsab/Prioritization/Plant_Cell/Cross_validation/10_fold/Results/Total_topology/we.RData")