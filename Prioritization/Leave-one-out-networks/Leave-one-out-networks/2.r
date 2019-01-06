load("~/prioritization/Combined_Network/New/Leave_one_out/Genemania.RData")
a1=a1_graph
library(Matrix)
library(igraph)
clos=closeness(a1,mode=c("all")) #Closeness 
save(clos,file="~/prioritization/Leave-one-out/New/Leave-one-out-network/2_Genemania.RData")

