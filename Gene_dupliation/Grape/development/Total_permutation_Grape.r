

library("isa2")
library("foreach")
library("doMC")
registerDoMC(20)
# expressionMat is the gene expression matrix

data = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/Vitis_vinifera/GE.mean.SD.txt", nrows = 23629, comment.char = "", header = T)
data = as.matrix(data)
# gene.names <- read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated CORNET2.0/Vitis_vinifera/gene.names.txt")
# gene.names = as.matrix(gene.names)

# Normalize gene expression in gene-wise and condition-wise
data.norm <- isa.normalize(data)
save(data.norm, file = "/scratch/iwt/ehsab/Gene_duplication/Results/SA/Local/Grape/data.norm.RData")

# For_loop for running 1000 permutations
foreach (i = 1:1000) %dopar%{
	print(i)
	sampled.Ec = t(apply(data.norm$Ec, 1, sample))
	sampled.Er = t(apply(data.norm$Er, 1, sample))
	save(sampled.Ec, sampled.Er, file = paste("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Local/Grape/Randomization/95/",i,".RData", sep = ""))
}


