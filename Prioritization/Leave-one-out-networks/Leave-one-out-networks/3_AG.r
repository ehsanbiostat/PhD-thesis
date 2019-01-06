load("~/prioritization/Combined_Network/New/Leave_one_out/AG.RData")
a1=a1_graph
library(Matrix)
library(igraph)
klein=unlist(authority.score(a1)[1]) #Klein authority index
save(klein,file="~/prioritization/Leave-one-out/New/Leave-one-out-network/3_AG.RData")

