

library(isa2)
library(foreach)
library(doMC)
registerDoMC(16)
load("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Input_data/Input_data.RData")

smartSeed <- matrix(0, nrow = dim(data.norm$Ec), ncol = 1) # every column is a seedvector

#
# Perform the ISA steps
#


filter <- function(x, t, dir) {
    if (dir=="updown") {
      x <- .Call("beta_filter_updown_vart", x, as.double(t), PACKAGE="isa2")
    } else if (dir=="up") {
      x <- .Call("beta_filter_up_vart", x, as.double(t), PACKAGE="isa2")
    } else { ## dir=="down"
      x <- .Call("beta_filter_down_vart", x, as.double(t), PACKAGE="isa2")
    }
  }


thr.row = 1.4 + (8 * 0.1)
gene_score_module = foreach (i = 1:dim(data.norm$Ec)[1], .combine = cbind) %dopar%{
	print(i)
	smartSeed[i, 1] <- 1
# 1) Calculate the condition scores for the seed gene ==> in case of only one seed gene this boils down to a vector of its gene expression values
	cond_scores <- data.norm$Er %*% smartSeed
	cond_scores <- apply(cond_scores, 2, function(x) x/sqrt(sum(x^2))) # normalize by the length of the vector ==> vector of length 1

# 2) Get the gene score vector
	gene_scores <- data.norm$Ec %*% cond_scores
	gene_scores <- apply(gene_scores, 2, function(x) x/sqrt(sum(x^2))) # normalize by the length of the vector ==> vector of length 1
	gene_scores <- filter(gene_scores, thr.row, "updown")

# 3) Recalculate the condition scores
	cond_scores <- data.norm$Er %*% gene_scores
	cond_scores <- apply(cond_scores, 2, function(x) x/sqrt(sum(x^2))) # normalize by the length of the vector ==> vector of length 1

	gene_score_module = ifelse(gene_scores > 0, 1, 0)
	# cond_score_module = ifelse(cond_scores > 0, 1, 0)
	
}


save(gene_score_module, file = paste("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Gene_score/Permutation/SA_threshold/Graph",thr.row,".RData", sep = ""))

