

## Reading Actual gene scores
list = read.table("~/Gene_duplication/Scripts/General/order.gene.score.txt")
list = as.vector(unlist(list))
setwd("/group/biocomp/users/ehsab/Gene_duplication/Results/SA/Gene_score/Total/")
SA_gene_score = lapply(list, read.table, sep= "\t" , header=T)


## Reading permuted gene scores
list = read.table("~/Gene_duplication/Scripts/General/order.thresh_gene_score.txt")
list = as.vector(unlist(list))
setwd("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Gene_score/Permutation/Threshold99/")
thresh_gene_score = lapply(list, read.table, sep= "\t")

gene_score_module = matrix(NA, length(SA_gene_score), length(SA_gene_score))
for (i in 1:length(SA_gene_score)){
	print(i)
	gene_score_module[, i] = ifelse(unlist(SA_gene_score[i]) > unlist(thresh_gene_score[i]), 1, 0)
}

