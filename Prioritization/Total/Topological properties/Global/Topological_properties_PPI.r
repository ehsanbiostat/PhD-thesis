a=read.table("~/prioritization/Combined_Network/New/Indivitual_Networks/PPI.txt")
a=a[,-3]
a=as.matrix(a)
a=rbind(a,matrix(c(27376,27376),1,2)) ##
  library(Matrix)
  a=sparseMatrix(a[,1],a[,2],x=1)
  a[dim(a)[1],dim(a)[1]]=0
  library(igraph)
  a=a+t(a)
  a1=graph.adjacency(a,mode=c("undirected"))
  graph.density(a1)
  transitivity(a1, type=c("global"))
  diameter(a1, directed = F, unconnected = T)
  average.path.length(a1, directed=F, unconnected=TRUE)