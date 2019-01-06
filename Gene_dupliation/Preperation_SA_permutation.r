

library("isa2")

# expressionMat is the gene expression matrix
data.SD = read.table("/group/biocomp/users/ehsab/Gene_duplication/Integrated CORNET2.0/Integrated.txt", colClasses = "numeric", nrows = 21428, comment.char = "", header = T)
gene.names <- read.table("/group/biocomp/users/ehsab/Gene_duplication/Integrated CORNET2.0/gene.names.txt")
gene.names = as.matrix(gene.names)
low.level.express = read.table("/group/biocomp/users/ehsab/Gene_duplication/Results/Expression_level/low.level.genes.SD.txt", header = F)

data.SD = as.matrix(data.SD)
## 19285 genes with 2933 conditions based on SD
data.SD = data.SD[-as.matrix(low.level.express), ]
gene.names = gene.names[-as.matrix(low.level.express)]


data.norm <- isa.normalize(data.SD)

#
# Prepare seed set 
#
## Reading 

save(gene.names, data.norm, file = "/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Input_data/Input_data.RData")
smartSeed <- matrix(0, nrow = nrow(data.SD), ncol = 1) # every column is a seedvector
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

write.table(gene_scores, file = "~/../../group/biocomp/users/ehsab/Gene_duplication/Results/SA/Gene_score/Total/G_Total.SD1.txt", row.names = F, col.names = T)
write.table(cond_scores, file = "~/../../group/biocomp/users/ehsab/Gene_duplication/Results/SA/Condition_score/Total/C_Total.SD1.txt", row.names = F, col.names = T)