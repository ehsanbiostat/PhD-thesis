load("~/prioritization/Combined_Network/New/Total/Total_list.RData")
list=read.table("~/prioritization/Combined_Network/Network_Codes.txt")
a=Total_list

i=4
# Remove ith netwrok from the "list"
a1=a[!a[,3]==list[i,2],]

# Sparse matrix
b=data.frame(27376,27376)
colnames(a1)=colnames(b)=c("A","B")
a1=rbind(a1[,-3],b)
library(Matrix)
a1_matrix=sparseMatrix(a1[,1],a1[,2],x=1)
a1_matrix[dim(a1_matrix)[1],dim(a1_matrix)[1]]=0

# graph object
library(igraph)
a1_matrix=a1_matrix+t(a1_matrix)
a1_graph=graph.adjacency(a1_matrix,mode=c("undirected"))

save(a1_graph,file="~/prioritization/Combined_Network/New/Leave_one_out/Mutant.RData")





