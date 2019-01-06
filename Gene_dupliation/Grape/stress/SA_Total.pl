
#!/usr/local/bin/perl -w

use strict;

for(my $i = 1; $i < 7180; $i++){
	
	my $jobfile = "SA".$i."_Tot.SD.r";
     open(OUT,">$jobfile");
     print OUT "

library(\"isa2\")

# expressionMat is the gene expression matrix
data = read.table(\"/ngsprojects/iwt/ehsab/Gene_duplication/Integrated CORNET2.0/Vitis_vinifera/Stress/Exp.matrix.txt\", nrows = 7179, comment.char = \"\", header = T)
gene.names <- rownames(data)
gene.names = as.matrix(gene.names)

data = as.matrix(data)


data.norm <- isa.normalize(data)

#
# Prepare seed set 
#
## Reading 


smartSeed <- matrix(0, nrow = nrow(data), ncol = 1) # every column is a seedvector
smartSeed[".$i.", 1] <- 1

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


colnames(cond_scores) = colnames(gene_scores) = gene.names[".$i."]

write.table(gene_scores, file = \"/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Compendium_wise/stress/Raw/Gene_score/G_Total.SD".$i.".txt\", row.names = F, col.names = T)
write.table(cond_scores, file = \"/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Compendium_wise/stress/Raw/Condition_score/C_Total.SD".$i.".txt\", row.names = F, col.names = T)"
	;
	close OUT;
}





