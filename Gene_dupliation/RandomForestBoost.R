

#################################################################
#
# Find other similar conditions to the top ten conditions for a given APs query
#
#################################################################

library(igraph)
library(foreach)
library(doMC)
library(randomForest)
options(width = 200)
library(pheatmap)
library(circlize)
# registerDoMC(2)

#################################################################################################
# Functions
#################################################################################################



Annot = function(x) {
	Annot = foreach(i = 1:length(x), .combine = c) %do% {
		strsplit(substring(x, 2), "[.]")[[i]][1]
	}
	return(Annot)
}


cores = 30
forkCL <- makeForkCluster(cores)
clusterSetRNGStream(forkCL, iseed = 123)
clusterEvalQ(forkCL, library(varSelRF))

RandomForestBoost = function(data, x, y, namex, namey, ntrees, classes, bootnumber, cores) {
	require(varSelRF)
	require(parallel)
	if(classes == 2) Data = rbind.data.frame(data[x,], data[y,])
	if(classes > 2) Data = data
	Class = as.factor(c(rep(namex, length(x)), rep(namey, length(y))))
	Result = varSelRFBoot(Data, Class, ntrees = ntrees, TheCluster = forkCL, bootnumber = bootnumber)
	return(Result)
	# stopCluster(forkCL)
}




#################################################################################################

k = 6
condition.names = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/condition_name.txt", header = T)
NAME = names(table(condition.names))
NAME.spe = NAME[k]
source("/ngsprojects/iwt/ehsab/Gene_duplication/Scripts/General/myImagePlot.R")
source("/ngsprojects/iwt/ehsab/Gene_duplication/Scripts/General/myImagePlot_sym.R")
ref = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/Reference_SD.txt", header = T)


# Expression data
data.exp = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/Integrated_SD.txt", colClasses = "numeric", nrows = 19285, comment.char = "", header = T)
data.exp.dev = data.exp[, condition.names == NAME.spe]
data.exp.scale = apply(data.exp.dev,1,scale)
data.exp.scale = t(data.exp.scale)
data.exp.scale = apply(data.exp.scale, 2, scale)
colnames(data.exp.scale) = colnames(data.exp.dev)
rownames(data.exp.scale) = ref[,1]
data.exp.scale.all = apply(data.exp,1,scale)
data.exp.scale.all = t(data.exp.scale.all)
data.exp.scale.all = apply(data.exp.scale.all, 2, scale)
colnames(data.exp.scale.all) = colnames(data.exp)
rownames(data.exp.scale.all) = ref[,1]
AA = Annot(colnames(data.exp.scale.all))
data.exp.scale.dev = data.exp.scale.all[, -which(is.na(as.numeric(AA)) | AA %in% Annot(colnames(data.exp.scale)) | condition.names == NAME.spe)]
data.exp.scale.dev = data.exp.scale.dev[,which(!duplicated(Annot(colnames(data.exp.scale.dev))))]


# Anchorpoints within two modules
AnchPollen = read.delim("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/Anchorpoint/Geneset1_2_53_90.txt.txt")
AnchRoot = read.delim("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/Anchorpoint/AT1G50090_AT3G19710_52_90.txt")
AnchRootPollen = read.delim("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/Anchorpoint/AnchRootPollen_90.txt")
AnchPollenRef = AnchPollen$Ref
AnchPollenLeafRef = AnchPollen$References
AnchRootRef = AnchRoot$Ref
AnchRootLeafRef = AnchRoot$References


#######################################################################################################################################################################
# Classification method to select Top conditions
#######################################################################################################################################################################



PollenLeaf = RandomForestBoost(classes = 2, data.exp.scale, as.character(AnchPollen$Pollen), as.character(AnchPollen$Leaf), "Pollen", "Leaf", ntrees = 20000, bootnumber = 1000, cores = cores)
RootLeaf = RandomForestBoost(classes = 2, data.exp.scale, as.character(AnchRoot$Root), as.character(AnchRoot$Leaf), "Root", "Leaf", ntrees = 20000, bootnumber = 1000, cores = cores)
RootPollen = RandomForestBoost(classes = 2, data.exp.scale, as.character(AnchRootPollen$Root.Pollen), as.character(AnchRootPollen$Leaf), "Root & Pollen", "Leaf", ntrees = 20000, bootnumber = 1000, cores = cores)

# PollenLeafAdd = RandomForestBoost(classes = 2, data.exp.scale.dev, as.character(AnchPollen$Pollen), as.character(AnchPollen$Leaf), "Pollen", "Leaf", ntrees = 20000, bootnumber = 1000, cores = cores)
# RootLeafAdd = RandomForestBoost(classes = 2, data.exp.scale.dev, as.character(AnchRoot$Root), as.character(AnchRoot$Leaf), "Root", "Leaf", ntrees = 20000, bootnumber = 1000, cores = cores)
# RootPollenAdd = RandomForestBoost(classes = 2, data.exp.scale.dev, as.character(AnchRootPollen$Root.Pollen), as.character(AnchRootPollen$Leaf), "Root & Pollen", "Leaf", ntrees = 20000, bootnumber = 1000, cores = cores)


#######################################################################################################################################################################

save(PollenLeaf, RootLeaf, RootPollen,
	file = "/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Plot/RandomForestBoost/All1.RData")





#######################################################################################################################################################################################
# For entire module (Pollen-Leaf & Root-Leaf)
#######################################################################################################################################################################################
Pollen = read.table("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/Module/Pollen_90.txt")
PollenLeaf = read.table("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/Module/Leaf_90.txt")
Root = read.table("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/Module/AT1G50090_Root_90.txt")
RootLeaf = read.table("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/Module/AT3G19710_Leaf_90.txt")


pollen = RandomForestBoost(classes = 2, data.exp.scale, as.character(Pollen[,1]), as.character(PollenLeaf[,1]), "Pollen", "Leaf", ntrees = 20000, bootnumber = 1000, cores = cores)
root = RandomForestBoost(classes = 2, data.exp.scale, as.character(Root[,1]), as.character(RootLeaf[,1]), "Root", "Leaf", ntrees = 20000, bootnumber = 1000, cores = cores)
#######################################################################################################################################################################################



save(pollen, root, PollenLeaf, RootLeaf, RootPollen, PollenLeafAdd, RootLeafAdd, RootPollenAdd,
	file = "/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Plot/RandomForestBoost/All.RData")
