load("~/prioritization/Combined_Network/New/Total/Total_list.RData")
c=read.table("~/prioritization/Combined_Network/New/Filtering_list/cregiyg.txt")
a=subset(Total_list,Total_list[,1]%in%c[,2])
a1=subset(Total_list,Total_list[,2]%in%c[,2])
a=rbind(a,a1)
a=a[!duplicated(a),]
library(igraph)
g=graph.data.frame(a[,-3],directed=F)



d=read.table("~/prioritization/Leave-one-out/New/Total_topology/global.txt",header=T)
d=d[c[,2],]
