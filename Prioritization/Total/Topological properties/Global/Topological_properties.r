load("~/prioritization/Combined_Network/New/Total/Total_graph.RData")
library(igraph)
a1=Total_graph
  #graph.density(a1)
  transitivity(a1, type=c("global"))
  #diameter(a1, directed = F, unconnected = T)
  average.path.length(a1, directed=F, unconnected=TRUE)