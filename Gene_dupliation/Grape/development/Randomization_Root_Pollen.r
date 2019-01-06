

	 
library(fields)
library(igraph)
library(foreach)
library(doMC)
registerDoMC(40)

condition.names = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/condition_name.txt", header = T)
NAME = names(table(condition.names))

k = 6
	
## Grape modules
load("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Graph_V2.RData")
Int = function(x,y){length(intersect(x, y))}	

PollenRootLeafclean = read.delim("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Orthology/Grape-Arabidopsis/PollenRootLeafclean.txt")
Triangle = read.delim("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Orthology/Grape-Arabidopsis/AllcleanZhen.txt")
gene.names = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/Vitis_vinifera/gene.names.txt")

SamplePollen = lapply(1:10000, function(x) {sample(as.character(Triangle$Name), length(unique(PollenRootLeafclean[PollenRootLeafclean$Module == "Pollen", "name"])))})
SampleRoot = lapply(1:10000, function(x) {sample(as.character(Triangle$Name), length(unique(PollenRootLeafclean[PollenRootLeafclean$Module == "Root", "name"])))})
	
# Score = matrix(NA, length(V(g)), 6)
# colnames(Score) = c("Grape", "Grape.size", "PollenLeaf", "PollenLeafPvalue", "RootLeaf", "RootLeafPvalue")
# Score = as.data.frame(Score)

# for (j in 1:length(V(g))){
# 	print(j)
# 	memb.comun.grape = data.frame("name" = gene.names[neighborhood(g, 1, j, mode = c("out"))[[1]],])
#     Rand = sapply(SamplePollen, function(x){sum(memb.comun.grape[, "name"] %in% x)})
#     PollenObse = sum(memb.comun.grape[, "name"] %in% PollenRootLeafclean[PollenRootLeafclean$Module == "Pollen", "Name"])
#     Score[j, "PollenLeafPvalue"] = sum(PollenObse <= Rand)/10000

# 	Rand = sapply(SampleRoot, function(x){sum(memb.comun.grape[, "name"] %in% x)})
#     RootObse = sum(memb.comun.grape[, "name"] %in% PollenRootLeafclean[PollenRootLeafclean$Module == "Root", "Name"])
#     Score[j, "RootLeafPvalue"] = sum(RootObse <= Rand)/10000

#     Score[j, "Grape"] = j
# 	Score[j, "Grape.size"] = nrow(memb.comun.grape)
# 	Score[j, "PollenLeaf"] = PollenObse
# 	Score[j, "RootLeaf"] = RootObse
# }


# Score = foreach (j = 1:length(V(g)), .combine = rbind) %dopar% {
# 	print(j)
# 	memb.comun.grape = data.frame("name" = gene.names[neighborhood(g, 1, j, mode = c("out"))[[1]],])
#     RandPollen = sapply(SamplePollen, function(x){sum(memb.comun.grape[, "name"] %in% x)})
#     PollenObse = sum(memb.comun.grape[, "name"] %in% PollenRootLeafclean[PollenRootLeafclean$Module == "Pollen", "name"])
    
# 	RandRoot = sapply(SampleRoot, function(x){sum(memb.comun.grape[, "name"] %in% x)})
#     RootObse = sum(memb.comun.grape[, "name"] %in% PollenRootLeafclean[PollenRootLeafclean$Module == "Root", "name"])
    

#     data.frame("Grape" = j, "Grape.size" = nrow(memb.comun.grape), "PollenLeaf" = PollenObse, "PollenLeafPvalue" = sum(PollenObse <= RandPollen)/10000,
# 		"RootLeaf" = RootObse, "RootLeafPvalue" = sum(RootObse <= RandRoot)/10000)
# }


Score = foreach (j = 1:length(V(g)), .combine = rbind) %dopar% {
	print(j)
	memb.comun.grape = data.frame("name" = gene.names[neighborhood(g, 1, j, mode = c("out"))[[1]],])
	memb.comun.grape.orth = sum(as.character(memb.comun.grape$name) %in% unique(as.character(Triangle$Name)))
	
	SampleGrape = lapply(1:10000, function(x) {sample(as.character(Triangle$Name), memb.comun.grape.orth)})
	RandPollen = mapply(Int, SamplePollen, SampleGrape)
    PollenObse = sum(memb.comun.grape[, "name"] %in% PollenRootLeafclean[PollenRootLeafclean$Module == "Pollen", "name"])
    
    SampleRootGrape = lapply(1:10000, function(x) {sample(as.character(Triangle$Name), memb.comun.grape.orth)})
	RandRoot = mapply(Int, SampleRoot, SampleGrape)
	RootObse = sum(memb.comun.grape[, "name"] %in% PollenRootLeafclean[PollenRootLeafclean$Module == "Root", "name"])
    

    data.frame("Grape" = j, "Grape.size" = nrow(memb.comun.grape), "ORTH.size" = memb.comun.grape.orth, "PollenLeaf" = PollenObse, "PollenLeafPvalue" = sum(PollenObse <= RandPollen)/10000,
		"RootLeaf" = RootObse, "RootLeafPvalue" = sum(RootObse <= RandRoot)/10000)
}





Score$PollenLeafFDR = p.adjust(Score$PollenLeafPvalue, method = "fdr")
Score$RootLeafFDR = p.adjust(Score$RootLeafPvalue, method = "fdr")

write.table(Score, file = "/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Compendium_wise/development/PollenRoot99.txt", row.names = F, quote = F, sep = "\t")





