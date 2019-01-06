
##------------------------------------------------------------------------------------------------------------------------------
# Reading Gene and condition scores from single file and integrate them to two matrices
##------------------------------------------------------------------------------------------------------------------------------
# change directory to location which GO.txt files are there
# Use ~/Gene_duplication/Scripts/General/order.pl to make a list for reading .txt files
list = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Scripts/General/order.Threshold.txt")
list = as.vector(unlist(list))

## Read text file with read.table function
setwd("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Compendium_wise/development/Score99/Threshold")
g = lapply(list, read.table, sep= "\t" , header=T)

cluster.size = dim(data.frame(g[1]))[1]
Threshold = matrix(NA, cluster.size * length(list), 2)
Threshold = data.frame("Threshold95" = NA, "Threshold99" = NA)

a = read.table("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Compendium_wise/development/Score99/Scores_V2/Score1.txt", header = T)
cluster.size = sort(unique(a$orthologous.cluster.size))[-c(1:3)]


for (i in 1:length(g)){
	print(i)
	Threshold = rbind(Threshold, data.frame(g[i]))
}
Threshold = Threshold[-1,]


Threshold = data.frame("Grape" = unlist(lapply(c(4:1607), rep, length(cluster.size))), "Arab" = rep(cluster.size, length(list)), Threshold)

save(Threshold, file = "/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Compendium_wise/development/Score99/Threshold/Threshold.RData")




































