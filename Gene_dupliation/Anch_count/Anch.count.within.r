

library(igraph)
library(foreach)
library(doMC)
library(data.table)
library(gtools)
registerDoMC(30)

## Read condition names
# condition.names = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/condition_name.txt", header = T)
# NAME = names(table(condition.names))
# k = 1
# NAME.spe = NAME[k]
# load(paste("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Gene_score/Compendium_wise/",NAME[k],"/Graph_V2.RData", sep = ""))
# paralogous = read.table("~/Gene_duplication/Seed_genes/paralogous_total.txt", header = T)
paralogous = read.table("~/Gene_duplication/Seed_genes/paralogous_NoFilter.txt", header = T)
# PARA = unique(paralogous$Ref)
# PARA.neigh = neighborhood(g, 1, PARA, mode = c("out"))
load("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Merge_module/abiotic/G1000.RData")
PARA.neigh = G2

## -----------------------------------------------------------------------------------------------------------------------------
# Count anchorpoint pairs within each module
## -----------------------------------------------------------------------------------------------------------------------------

# NO.anch.within = foreach (j = 1:length(PARA.neigh), .combine = rbind.data.frame) %dopar% {
# 	print(j)
	
# 	RR = foreach (i = 1:(nrow(paralogous)), .combine = c) %do% {
# 		sum(unlist(PARA.neigh[j]) %in% as.numeric(paralogous[i,1:2])) == 2
# 	}

# 	data.frame("PARA" = j, "Count" = sum(RR)/2, "Recent" = sum(paralogous[RR, "WGDevent"] == "recent", na.rm = T)/2, 
# 		"Ancient" = sum(paralogous[RR, "WGDevent"] == "ancient", na.rm = T)/2, "module.size" = length(unlist(PARA.neigh[j])))
# }

# NO.anch.within["Anch.perc"] = NO.anch.within$Count/NO.anch.within$module.size
# # NO.anch.within = NO.anch.within[order(-NO.anch.within$Anch.perc), ]
# write.table(NO.anch.within, file = "/scratch/iwt/ehsab/Gene_duplication/Results/SA/Anch_count/Anch.count.within_1000.txt", row.names = F, sep = "\t", quote = F)

## -----------------------------------------------------------------------------------------------------------------------------



## -----------------------------------------------------------------------------------------------------------------------------
# Permutation test
## -----------------------------------------------------------------------------------------------------------------------------
# NO.anch.within = read.delim("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Anch_count/Anch.count.within.txt")
# modu.size = sort(unique(NO.anch.within$module.size))[-1]
# Permut = matrix(NA, 10000, length(modu.size))

# # To use values of ‘n’ above about 45, you will need to increase R's
# # recursion limit.  See the ‘expression’ argument to the ‘options’
# # command for details on how to do this.
# options(expressions=1e5)

# for (j in 1:length(modu.size)) {
# 	print(j)
# 	Perm.anch.within = foreach (k = 1:10000, .combine = c) %dopar% {
# 		x = sample(1:19285, modu.size[j])
# 		# RR = foreach (i = 1:(nrow(paralogous)), .combine = c) %do% {
# 		# 	sum(x %in% as.numeric(paralogous[i,1:2])) == 2
# 		# }
# 		# sum(RR)/2
# 		# Generate all combination of randomly sampled data as a two columns matrix
# 		# Then comparing the random generated list with paralogous list is 50 time faster than the
# 		# conventional way
# 		M1 = setkey(data.table(combinations(length(x), 2, v = x)))
# 		M2 = setkey(data.table(paralogous[, 1:2]))
# 		length(na.omit(M2[M1,which=TRUE]))
# 	}
# 	Permut[,j] = Perm.anch.within
# }

# save(Permut, file = "/scratch/iwt/ehsab/Gene_duplication/Results/SA/Anch_count/Permut.RData")

# apply(Permut, 2, quantile, probs = 0.95)




## -----------------------------------------------------------------------------------------------------------------------------
# Count co-appearance of each anchorpoint pairs within all modules
## -----------------------------------------------------------------------------------------------------------------------------

coappearance.anch.within = foreach (i = 1:nrow(paralogous), .combine = rbind.data.frame) %dopar% {
	print(i)
	
	RR = foreach (j = 1:length(PARA.neigh), .combine = c) %do% {
		sum(unlist(PARA.neigh[j]) %in% as.numeric(paralogous[i,1:2]))
	}

	data.frame(paralogous[i,], "Count" = sum(RR == 1), "co.appearance" = sum(RR == 2))
}

write.table(coappearance.anch.within, file = "/scratch/iwt/ehsab/Gene_duplication/Results/SA/Anch_count/coappearance.anch.within_1000.txt", row.names = F, sep = "\t", quote = F)
coappearance.anch.within["Rate"] = coappearance.anch.within$co.appearance / coappearance.anch.within$Count
coappearance.anch.within.zero = coappearance.anch.within[coappearance.anch.within$co.appearance == 0,]
coappearance.anch.within.zero = coappearance.anch.within.zero[order(-coappearance.anch.within.zero$Count), ]
coappearance.anch.within = coappearance.anch.within[order(-coappearance.anch.within$Rate, -coappearance.anch.within$co.appearance), ]
coappearance.anch.within.inf = coappearance.anch.within[coappearance.anch.within$Rate != Inf, ]
## -----------------------------------------------------------------------------------------------------------------------------

memb.comun = unique(coappearance.anch.within.zero[1:200, "Ref"])









