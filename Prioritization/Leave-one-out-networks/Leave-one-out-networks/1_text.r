load("~/prioritization/Combined_Network/New/Leave_one_out/text.RData")
a1=a1_graph
library(Matrix)
library(igraph)
alldeg=degree(a1)        #Total number of connections
bet=betweenness(a1,directed=F)  # Betweennees Centrality
save(alldeg,bet,file="~/prioritization/Leave-one-out/New/Leave-one-out-network/1_text.RData")

