load("/ngsprojects/iwt/ehsab/Prioritization/prioritization/Combined_Network/New/Total/Total_matrix.RData")
load("/ngsprojects/iwt/ehsab/Prioritization/prioritization/Combined_Network/New/Total/Total_graph.RData")
cregiyg = read.table("/ngsprojects/iwt/ehsab/Prioritization/Plant_Cell/Cross_validation/10_fold/Seed_genes/10.txt")
cregiyg = cregiyg[, 2]
a = Total_matrix
a1 = Total_graph
library(Matrix)
library(igraph)
p = a[, unlist(cregiyg)]
d = p %*% t(p)   #Matrix multiplication shared genes with IYG
sh = d[, unlist(cregiyg)]
IYG_S = apply(sh, 1, sum)
ciyg=a[, unlist(cregiyg)]
decregiyg=apply(ciyg, 1, sum) #Number of connection with IYG & REG
# alldeg=degree(a1)        #Total number of connections
cregiyg=as.matrix(cregiyg)
# bet=betweenness(a1,directed=F)  # Betweennees Centrality
# trans=transitivity(a1,type=c("local"),isolates=c("zero")) #Clutering Coefficient
# clos=closeness(a1,mode=c("all")) #Closeness 
# klein=unlist(authority.score(a1)[1]) #Klein authority index

sim.jacard=similarity.jaccard(a1) #
sim.jacard=sim.jacard[,unlist(cregiyg)]
sim.jacard=apply(sim.jacard,1,sum)

sim.dic=similarity.dice(a1)
sim.dic=sim.dic[,unlist(cregiyg)]
sim.dic=apply(sim.dic,1,sum)


sim.we=similarity.invlogweighted(a1)
sim.we=sim.we[,unlist(cregiyg)]
sim.we=apply(sim.we,1,sum)


short=apply(shortest.paths(a1,c(1:dim(a)[1]),c(unlist(cregiyg))),1,sum)


save(alldeg,deiyg,dereg,IYG_S,REG_S,bet,trans,clos,klein,sim.jacard_iyg,sim.jacard_reg,sim.dic_iyg,sim.dic_reg,sim.we_iyg,sim.we_reg,short_iyg,short_reg,file="ren1.RData")