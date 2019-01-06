
library(igraph)
library(foreach)
library(doMC)
registerDoMC(3)

## Read condition names
condition.names = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/condition_name.txt", header = T)
NAME = names(table(condition.names))
k = 1
NAME.spe = NAME[k]

## Read text file with read.table function
setwd(paste("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Local/Compendium_wise/",NAME[k],"/Score99/Anch/Scores_V2/", sep = ""))

list = list.files()
list = as.vector(unlist(list))

g = lapply(list, read.table, sep= "\t" , header=T)
Total = do.call(rbind.data.frame, g)

Total.25 = Total[Total$Jaccard > 0.24,]
Total.25 = Total.25[!duplicated(Total.25[, c("Paralogous.cluster", "Cluster.overlap", "Normalized.Paralogous.cluster")]),]
Total.25 = Total.25[order(-Total.25$Jaccard), ]
Total.25 = Total.25[!mapply(function(x,y){x == y}, x = Total.25$Gene, y = Total.25$Paralogous),]
# paralogous = read.table("~/Gene_duplication/Seed_genes/paralogous_total.txt", header = T)
paralogous = read.table("~/Gene_duplication/Seed_genes/paralogous_NoFilter.txt", header = T)
PARA = unique(paralogous$Ref)


# --------------------------------------------------------------------------------------------------------------------------------------------------
## Condition scores
# --------------------------------------------------------------------------------------------------------------------------------------------------
list = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Scripts/General/order.condition_score.txt")
list = as.vector(unlist(list))
## Read text or RData file with read.table or load function
setwd(paste("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Condition_score/Compendium_wise/",NAME.spe,"/Raw", sep = ""))
# C = lapply(list, read.table, sep= "\t" , header=T)
C = lapply(list, function(x) mget(load(x)))
total_cond_scores = matrix(NA, length(unlist(C[1])), 19285)

for (i in 1:19285){
	print(i)
	total_cond_scores[, i] = as.matrix(unlist(C[i]))
}
# --------------------------------------------------------------------------------------------------------------------------------------------------




##------------------------------------------------------------------------------------------------------------------------------
# Reading Gene scores from single file and integrate them to a matrix
##------------------------------------------------------------------------------------------------------------------------------
# change directory to location which Gene scores files are there
setwd(paste("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Gene_score/Compendium_wise/",NAME.spe,"/Raw", sep = ""))

## Generate name of Gene score files by a Perl script and read it to an object
list = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Scripts/General/order.gene.score.txt")
list = as.vector(unlist(list)[1:19285])

## Read text file with read.table function
g = lapply(list, function(x) mget(load(x)))

total_gene_scores = matrix(NA, length(unlist(g[1])), length(g))

for (i in 1:length(g)){
	print(i)
	total_gene_scores[, i] =  as.matrix(unlist(g[i]))
}
##------------------------------------------------------------------------------------------------------------------------------




# --------------------------------------------------------------------------------------------------------------------------------------------------
## Gene expression
data.exp = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/Integrated_SD.txt", colClasses = "numeric", nrows = 19285, comment.char = "", header = T)
data.exp = data.exp[, condition.names == NAME.spe]
load(paste("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Gene_score/Compendium_wise/",NAME[k],"/Graph_V2.RData", sep = ""))
condition_Jaccard = foreach (j = 1:nrow(Total.25), .combine = c) %dopar% {	
	# m = m + 20
 	print(j)
	memb.comun = neighborhood(g, 1, Total.25$Gene[j], mode = c("out"))[[1]]
	memb.comun.para = neighborhood(g, 1, Total.25$Paralogous[j], mode = c("out"))[[1]]
	
	CondRank1 = foreach(i = 1:length(memb.comun), .combine = rbind) %do%{
		rank(total_cond_scores[, memb.comun[i]])
	}

	CondRank2 = foreach(i = 1:length(memb.comun.para), .combine = rbind) %do%{
		rank(total_cond_scores[, memb.comun.para[i]])
	}

	
	CondRank = rbind(c(1:dim(total_cond_scores)[1]), CondRank1)
	colnames(CondRank) = colnames(data.exp)
	ORDER = 1
	CondRank0 = CondRank[,order(CondRank[ORDER+1, ], decreasing = T)]
	
	CondRank = rbind(c(1:dim(total_cond_scores)[1]), CondRank2)
	colnames(CondRank) = colnames(data.exp)
	ORDER = 1
	CondRank01 = CondRank[,order(CondRank[ORDER+1, ], decreasing = T)]
	length(intersect(colnames(CondRank0)[1:10], colnames(CondRank01)[1:10]))/length(union(colnames(CondRank0)[1:10], colnames(CondRank01)[1:10]))
	
}


Total.25["condition_Jaccard"] = condition_Jaccard
Total.25 = Total.25[order(-Total.25$Jaccard, -Total.25$condition_Jaccard), ]
Total.25 = Total.25[Total.25$condition_Jaccard > 0.5, ]
save(Total.25, Total, file = "/scratch/iwt/ehsab/Gene_duplication/Results/SA/Merge_module/abiotic/G1.RData")
PARA = unique(Total.25$Ref)

Dist = as.matrix(dist(t(total_cond_scores[,PARA])))
colnames(Dist) = rownames(Dist) = PARA

condition_dist = foreach (i = 1:nrow(Total.25), .combine = c) %dopar% {
	Dist[rownames(Dist) %in% Total.25$Gene[i], colnames(Dist) %in% Total.25$Paralogous[i]]
}

Total.25["condition_dist"] = condition_dist
Total.25["dist"] = condition_dist



# --------------------------------------------------------------------------------------------------------------------------------------------------
# Merging gene and condition score for all genes
# We need to read condition and gene score matrices, then merge them and apply a distance function
# --------------------------------------------------------------------------------------------------------------------------------------------------

total = cbind(total_gene_scores, t(total_cond_scores))
Dist = dist(total[PARA,])
Hclust = hclust(Dist, method = "ave")


# Total number of merged modules
Module.NO = 1000
G2.group = data.frame("seed" = PARA, "module" = cutree(Hclust, k = Module.NO))

load("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Merge_module/abiotic/G1.RData")
load(paste("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Gene_score/Compendium_wise/",NAME[k],"/Graph_V2.RData", sep = ""))
# Merge modules
G2 = foreach (i = 1:Module.NO) %dopar% {
	unique(unlist(neighborhood(g, 1, G2.group[G2.group$module == i, "seed"], mode = c("out"))))
}

names(G2) = sapply(1:Module.NO, function(i){paste(G2.group[G2.group$module == i, "seed"], collapse = "_")})

# names(G2) = mapply(function(x,y){paste(x,y,sep = "_")}, x = Total.25[,1], y = Total.25[,2])



## -----------------------------------------------------------------------------------------------------------------------------
# Count shared genes (overlap) and shared anchorpoints
## -----------------------------------------------------------------------------------------------------------------------------
G2.ov = foreach(i = 1:length(G2), .combine = rbind.data.frame) %dopar% {
	print(i)
	data.frame("Overlap" = mapply(function(x,y){length(intersect(x,y))}, x = G2, y = G2[i]), 
		"Para.overl" = mapply(function(x,y){length(intersect(paralogous[paralogous$Ref %in% x, "Anch"],y))}, x = G2, y = G2[i]),
		"NameA" = names(G2)[i], 
		"sizeB" = sapply(G2, length), "sizeA" = length(unlist(G2[i])))
}
rownames(G2.ov) = 1:nrow(G2.ov)
G2.ov["Jaccard"] = G2.ov$Overlap/(G2.ov$sizeA + G2.ov$sizeB)
G2.ov = G2.ov[order(-G2.ov$Jaccard),]
G2.ov = G2.ov[-c(1:Module.NO),]
G2.ov = G2.ov[order(-G2.ov$Para.overl),]
save(G2, G2.ov , Hclust, Total, file = "/scratch/iwt/ehsab/Gene_duplication/Results/SA/Merge_module/abiotic/G2.RData")
## -----------------------------------------------------------------------------------------------------------------------------




## -----------------------------------------------------------------------------------------------------------------------------
# Count anchorpoint pairs within each module
## -----------------------------------------------------------------------------------------------------------------------------


NO.anch.within = foreach (j = 1:length(PARA.neigh), .combine = rbind.data.frame) %dopar% {
	print(j)
	
	RR = foreach (i = 1:(nrow(paralogous)/2), .combine = c) %do% {
		sum(unlist(PARA.neigh[j]) %in% as.numeric(paralogous[i,1:2])) == 2
	}

	data.frame("Count" = sum(RR), "Recent" = sum(paralogous[RR, "WGDevent"] == "recent", na.rm = T)/2, 
		"Ancient" = sum(paralogous[RR, "WGDevent"] == "ancient", na.rm = T)/2, "module.size" = length(unlist(PARA.neigh[j])))
}

NO.anch.within["Anch.perc"] = NO.anch.within$Count/NO.anch.within$module.size
NO.anch.within[""] = 
NO.anch.within[""] = 

NO.anch.within = NO.anch.within[order(-NO.anch.within$Anch.perc), ]
## -----------------------------------------------------------------------------------------------------------------------------




## -----------------------------------------------------------------------------------------------------------------------------
# Fish out significant modules based on permuation test
## logic variable shows which GO.ov rows meet the permutation test conditions to be not-random
## -----------------------------------------------------------------------------------------------------------------------------
load("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Merge_module/abiotic/G2.RData")
logic = foreach (i = 1:nrow(G2.ov), .combine = c) %dopar% {
	k = which(Perm.test$sizeA %in% G2.ov$sizeA[i]  & Perm.test$sizeB %in% G2.ov$sizeB[i])
	G2.ov$Overlap[i] < Perm.test[k, "Threshold_Jac95"] & G2.ov$Para.overl[i] > Perm.test[k, "Threshold95"]
}
G2.ov[logic,]
## -----------------------------------------------------------------------------------------------------------------------------





# ----------------------------------------------------------------------------------------------------------------------------------------------- 
# Permutation test
# ----------------------------------------------------------------------------------------------------------------------------------------------- 
cluster.size.one = sort(unique(c(G2.ov$sizeA, G2.ov$sizeB)))
Perm.test = foreach (j = 1:length(cluster.size.one), .combine=rbind) %dopar%{
	print(j)
	cluster.size.Rand = rep(cluster.size.one[j], 1000)
	# Sample from nodes with the cluster's size
	cluster = lapply(cluster.size.Rand, function(x){sample(1:19285, x)})
	Anch.one = lapply(cluster, function(x) {paralogous[paralogous$Ref %in% x, "Anch"]})
	Int = function(x,y){length(intersect(x, y))}
	Result = foreach (i = 1:length(cluster.size.one), .combine=rbind) %do%{
		cluster.size.two = rep(cluster.size.one[i], 1000)
		Anch.two = lapply(cluster.size.two, function(x){sample(1:19285, x)})							
		Anc = mapply(Int, Anch.one, Anch.two)
		Jac = mapply(Int, cluster, Anch.two)
		data.frame("Threshold95" = quantile(Anc, probs = 0.95), "Threshold99" = quantile(Anc, probs = 0.99), 
			"Threshold_Jac95" = quantile(Jac, probs = 0.05), "Threshold_Jac99" = quantile(Jac, probs = 0.01), "sizeA" = cluster.size.one[i],
			 "sizeB" = cluster.size.one[j])
	}
	Result
}
# ----------------------------------------------------------------------------------------------------------------------------------------------- 




library(NbClust)
res = NbClust(x, diss=diss_matrix, distance = NULL, min.nc=2, max.nc=6,
                 method = "ward.D", index = "ch")








