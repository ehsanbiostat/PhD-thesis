AG_list=read.table("~/../../group/biocomp/users/ehsab/Combined_Network/Indivitual_Networks/AGRIS.txt")
Genemania_list=read.table("~/../../group/biocomp/users/ehsab/Combined_Network/Indivitual_Networks/Genemania.txt")
Genie3_list=read.table("~/../../group/biocomp/users/ehsab/Combined_Network/Indivitual_Networks/Genie3.txt")
PPI_list=read.table("~/../../group/biocomp/users/ehsab/Combined_Network/Indivitual_Networks/PPI.txt")
PCC_list=read.table("~/../../group/biocomp/users/ehsab/Combined_Network/Indivitual_Networks/PCC.txt")
Mutant_list=read.table("~/../../group/biocomp/users/ehsab/Combined_Network/Indivitual_Networks/Mutant.txt")
text_list=read.table("~/../../group/biocomp/users/ehsab/Combined_Network/Indivitual_Networks/text.txt")

# 27376 nodes - 24361406 edges
# list
Total_list=rbind(AG_list[,c(1,2)],Genemania_list[,c(1,2)],Genie3_list[,c(1,2)],Mutant_list[,c(1,2)],PCC_list[,c(1,2)],PPI_list[,c(1,2)],text_list[,c(1,2)])
list=c(rep(1,dim(AG_list)[1]),rep(2,dim(Genemania_list)[1]),rep(3,dim(Genie3_list)[1]),rep(4,dim(Mutant_list)[1]),rep(5,dim(PCC_list)[1]),rep(6,dim(PPI_list)[1]),rep(7,dim(text_list)[1]))
Total_list=cbind(Total_list,list)
save(Total_list,file="~/../../group/biocomp/users/ehsab/Combined_Network/Total_list.RData")

# Sparse matrix
a=data.frame(27376,27376)
colnames(Total_list)=colnames(a)=c("A","B")
Total_list=rbind(Total_list[,-3],a)
library(Matrix)
Total_matrix=sparseMatrix(Total_list[,1],Total_list[,2],x=1)
Total_matrix[dim(Total_matrix)[1],dim(Total_matrix)[1]]=0
save(Total_matrix,file="~/../../group/biocomp/users/ehsab/Combined_Network/Total_matrix.RData")


# graph object
library(igraph)
Total_matrix=Total_matrix+t(Total_matrix)
Total_graph=graph.adjacency(Total_matrix,mode=c("undirected"))
save(Total_graph,file="~/../../group/biocomp/users/ehsab/Combined_Network/Total_graph.RData")