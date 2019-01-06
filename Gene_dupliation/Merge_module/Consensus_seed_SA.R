

options(width = 200)
library(foreach)
library(doMC)
registerDoMC(3)

paralogous = read.delim("~/Gene_duplication/Seed_genes/paralogous_NoFilter.txt")
PARA = unique(paralogous$Ref)

# Concatenate all consensus cluster membership from all step 2 to 50 cluster
M = foreach (j = 2:100, .combine = rbind) %dopar% {
	print(j)
	a100 = read.csv(paste("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Merge_module/abiotic/untitled_consensus_cluster/untitled_consensus_cluster.k=",j,".consensusClass.csv", sep = ""), header = F)
	# a100 = read.csv(paste("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Merge_module/abiotic/200_5/200_5.k=",j,".consensusClass.csv", sep = ""), header = F)
	a100[,2]
}

# Count how many time a gene cluster membership has changed through cluster sizes froom 2 to 50
RR = lapply(1:ncol(M), function(x){table(M[,x])})
# Select the cluster which a gene has most frequently has been assigned to
RR.max = sapply(RR, function(x){names(x)[which(x == max(x))]})
# Select genes which has been assigned constantly more than 40 times to the same cluster through cluster sizes from 2 to 50 (80%)
sapply(RR, max) > 39

# Select stable paralogous genes and their corresponding clusters as seed genes for next step
G = data.frame("PARA" = PARA[which(sapply(RR, max) > 39)], "Cluster" = as.numeric(unlist(RR.max[sapply(RR, max) > 39])))
write.table(G, file = "/scratch/iwt/ehsab/Gene_duplication/Results/SA/Merge_module/abiotic/Seed.txt", row.names = F, sep = "\t")


