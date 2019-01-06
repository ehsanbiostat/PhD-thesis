
#!/usr/local/bin/perl -w

use strict;

for(my $i = 1; $i < 6588; $i++){
	
	my $jobfile = "SA".$i.".r";
     open(OUT,">$jobfile");
     print OUT "

library(\"isa2\")

# expressionMat is the gene expression matrix
data = read.table(\"~/../../group/biocomp/users/ehsab/Gene_duplication/Integrated CORNET2.0/Integrated.txt\", colClasses = \"numeric\", nrows = 21428, comment.char = \"\", header = T)
gene.names <- read.table(\"~/../../group/biocomp/users/ehsab/Gene_duplication/Integrated CORNET2.0/gene.names.txt\")
gene.names = as.matrix(gene.names)


data = as.matrix(data)
data.norm <- isa.normalize(data)

#
# Prepare seed set 
#
## Reading 
anchorpoint = read.table(\"~/Gene_duplication/Seed_genes/anchorpointTable_final.txt\", header = T)

## Combining homologous genes to a vector 
seedList = c(as.matrix(anchorpoint\$gene_x), as.matrix(anchorpoint\$gene_y))
## 7484 unique seed genes
seedList = sort(unique(seedList))


#seedList <- toupper(seedList)
## 6587 of seed genes are available in CORNET2.0 data set
seedIdx <- which(gene.names %in% seedList)

smartSeed <- matrix(0, nrow = nrow(data), ncol = 1) # every column is a seedvector
smartSeed[seedIdx[".$i."], 1] <- 1

#
# Perform the ISA steps
#

# 1) Calculate the condition scores for the seed gene ==> in case of only one seed gene this boils down to a vector of its gene expression values
	cond_scores <- data.norm\$Er %*% smartSeed
	cond_scores <- apply(cond_scores, 2, function(x) x/sqrt(sum(x^2))) # normalize by the length of the vector ==> vector of length 1

# 2) Get the gene score vector
	gene_scores <- data.norm\$Ec %*% cond_scores
	gene_scores <- apply(gene_scores, 2, function(x) x/sqrt(sum(x^2))) # normalize by the length of the vector ==> vector of length 1

# 3) Recalculate the condition scores
	cond_scores <- data.norm\$Er %*% gene_scores
	cond_scores <- apply(cond_scores, 2, function(x) x/sqrt(sum(x^2))) # normalize by the length of the vector ==> vector of length 1


colnames(cond_scores) = colnames(gene_scores) = seedList[".$i."]

write.table(gene_scores, file = \"~/../../group/biocomp/users/ehsab/Gene_duplication/Results/SA/Gene_score/G".$i.".txt\", row.names = F, col.names = T)
write.table(cond_scores, file = \"~/../../group/biocomp/users/ehsab/Gene_duplication/Results/SA/Condition_score/C".$i.".txt\", row.names = F, col.names = T)"
	;
	close OUT;
}





