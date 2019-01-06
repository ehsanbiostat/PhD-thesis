

# change directory to location which GO.txt files are there
# Use prioritization/New_IYG_Results/order.pl to make a list for reading GO.txt files
list=read.table("~/prioritization/New_IYG_Results/order.txt")
list=as.vector(unlist(list))
g=lapply(list,read.table,sep="\t",header=T)
save(g, file="GO.RData")



