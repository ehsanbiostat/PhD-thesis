

library("isa2")
library("foreach")
library("doMC")
registerDoMC(2)
# expressionMat is the gene expression matrix
data.exp = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated CORNET2.0/Integrated_SD.txt", colClasses = "numeric", nrows = 19285, comment.char = "", header = T)
condition.names = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated CORNET2.0/condition_name.txt", header = T)
NAME = names(table(condition.names))

gene.names <- read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated CORNET2.0/gene.names.low.mod.SD.txt")
gene.names = as.matrix(gene.names)

foreach (j = 1:length(NAME)) %dopar% {
	data = as.matrix(data.exp[, condition.names == NAME[j]])
	data.norm <- isa.normalize(data)

	#
	# Prepare seed set 
	#
	## Reading 


	smartSeed <- matrix(0, nrow = nrow(data), ncol = 1) # every column is a seedvector
	smartSeed[1, 1] <- 1



	#
	# Perform the ISA steps
	#

	# 1) Calculate the condition scores for the seed gene ==> in case of only one seed gene this boils down to a vector of its gene expression values
		cond_scores <- data.norm$Er %*% smartSeed
		cond_scores <- apply(cond_scores, 2, function(x) x/sqrt(sum(x^2))) # normalize by the length of the vector ==> vector of length 1

	# 2) Get the gene score vector
		gene_scores <- data.norm$Ec %*% cond_scores
		gene_scores <- apply(gene_scores, 2, function(x) x/sqrt(sum(x^2))) # normalize by the length of the vector ==> vector of length 1

	# 3) Recalculate the condition scores
		cond_scores <- data.norm$Er %*% gene_scores
		cond_scores <- apply(cond_scores, 2, function(x) x/sqrt(sum(x^2))) # normalize by the length of the vector ==> vector of length 1


	colnames(cond_scores) = colnames(gene_scores) = gene.names[1]

	# write.table(gene_scores, file = paste("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Gene_score/Compendium_wise/",as.name(NAME[j]),"/Raw/G1.txt", sep = ""), row.names = F, col.names = T)
	# write.table(cond_scores, file = paste("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Condition_score/Compendium_wise/",as.name(NAME[j]),"/Raw/C1.txt", sep = ""), row.names = F, col.names = T)

	save(gene_scores, file = paste("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Gene_score/Compendium_wise/",as.name(NAME[j]),"/Raw/G1.RData", sep = ""))
	save(cond_scores, file = paste("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Condition_score/Compendium_wise/",as.name(NAME[j]),"/Raw/C1.RData", sep = ""))
	
	}
	
	