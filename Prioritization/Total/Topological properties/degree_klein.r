load("/ngsprojects/iwt/ehsab/Prioritization/prioritization/Combined_Network/New/Total/Total_matrix.RData")
load("/ngsprojects/iwt/ehsab/Prioritization/prioritization/Combined_Network/New/Total/Total_graph.RData")
cregiyg = read.table("/ngsprojects/iwt/ehsab/Prioritization/Plant_Cell/Cross_validation/10_fold/Seed_genes/10.txt")
cregiyg=cregiyg[, 2]
a = Total_matrix
a1 = Total_graph
library(Matrix)
library(igraph)

p= a[, unlist(cregiyg)]
d = p %*% t(p)   #Matrix multiplication shared genes with IYG
sh = d[, unlist(cregiyg)]
cregiyg_sh = apply(sh ,1 ,sum)
ciyg = a[, unlist(cregiyg)]
decregiyg = apply(ciyg ,1 ,sum) #Number of connection with IYG & REG
alldeg=degree(a1)        #Total number of connections
cregiyg = as.matrix(cregiyg)
klein=unlist(authority.score(a1)[1]) #Klein authority index

save(alldeg,decregiyg,cregiyg_sh,klein, file = "/ngsprojects/iwt/ehsab/Prioritization/Plant_Cell/Cross_validation/10_fold/Results/Total_topology/degree_klein.RData")