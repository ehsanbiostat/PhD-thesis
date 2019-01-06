

library("isa2")
library("foreach")
library("doMC")
registerDoMC(30)
# expressionMat is the gene expression matrix

data.exp = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/Integrated_SD.txt", colClasses = "numeric", nrows = 19285, comment.char = "", header = T)
condition.names = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/condition_name.txt", header = T)
NAME = names(table(condition.names))

gene.names <- read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/gene.names.low.mod.SD.txt")
gene.names = as.matrix(gene.names)


for (j in 1:1) {
	
	print(paste(as.name(NAME[j])))
	data = as.matrix(data.exp[, condition.names == NAME[j]])
	# Normalize gene expression in gene-wise and condition-wise
	data.norm <- isa.normalize(data)

	# save(data.norm, file = paste("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Local/Compendium_wise/",as.name(NAME[j]),"/data.norm.RData", sep = ""))

	# For_loop for running 1000 permutations
	foreach (i = 1:1000) %dopar%{
		print(i)

		sampled.Ec = t(apply(data.norm$Ec, 1, sample))
		sampled.Er = t(apply(data.norm$Er, 1, sample))
		
		save(sampled.Ec, sampled.Er, file = paste("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Local/Compendium_wise/",as.name(NAME[j]),"/Randomization/95/",i,".RData", sep = ""))
	}
}

