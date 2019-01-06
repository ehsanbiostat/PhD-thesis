
##########################################################################
# Libraries
##########################################################################

library("isa2")
library("foreach")
library("doMC")
library(igraph)
library(pheatmap)
library(preprocessCore)
options(width = 200)
source("/ngsprojects/iwt/ehsab/Gene_duplication/Scripts/General/myImagePlot.R")
source("/ngsprojects/iwt/ehsab/Gene_duplication/Scripts/General/myImagePlot_sym.R")


##########################################################################
# Functions
##########################################################################

# Scaling a matrix on the row and column
source("/ngsprojects/iwt/ehsab/Gene_duplication/Scripts/Functions/Scaling.R")
# Heatmaps by pheatmap package
source("/ngsprojects/iwt/ehsab/Gene_duplication/Scripts/Functions/HeatMap.R")
# Tissue-specific score calculation 
source("/ngsprojects/iwt/ehsab/Gene_duplication/Scripts/Functions/EEI.R")
# Correlation between two rows of two different matrices
source("/ngsprojects/iwt/ehsab/Gene_duplication/Scripts/Functions/CORR.R")

# Read Marcopolo output data
readMarco = function(marco, Module) {
	colnames(marco) = c(as.matrix(marco[2,]))
	MM = marco[-c(1:5), -c(1:3)]
	Data = t(apply(MM, 1, as.numeric))
	colnames(Data) = colnames(MM)
	rownames(Data) = marco[-c(1:5), 1]
	Data = cbind.data.frame(Data, "Unique.ID" =  rownames(Data))
	rownames(Data) = as.character(merge(Data, map, by = "Unique.ID", sort = F)$name)
	Data = Data[, -ncol(Data)]
	if(Module != "conserve") Data = cbind.data.frame(Data, "class" = Module)
	if(Module == "conserve") {
		PolleN = read.delim("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Grape/Multiple_othologous_modules/Pollen.txt")
		RooT = read.delim("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Grape/Multiple_othologous_modules/Root.txt")
		Lable = rbind(PolleN[,c("Grape", "Module")], RooT[,c("Grape", "Module")])
		Data = cbind.data.frame(Data, "Grape" =  rownames(Data))
		Data = merge(Data, Lable, by = "Grape")
		Data = Data[!duplicated(Data[,1]),]
		rownames(Data) = Data[,1]
		Data = Data[,-1]
		colnames(Data)[ncol(Data)] = "class"
	}
	return(Data)
}


# Draw heatmap with ordered conditions based on median condition scores of genes within the module (j) with image score
OrderCondition = function(j, heatmap, Module, clusterRow, clusterCol, cutrow, cutcol, Top, address, binary) {
	CondRank = foreach(i = 1:length(j), .combine = rbind) %do%{
		rank(total_cond_scores[, j[i]])
	}
	Med.cond.score = apply(CondRank, 2, median)
	CondRank = rbind(c(1:dim(total_cond_scores)[1]), Med.cond.score, CondRank)
	colnames(CondRank) = colnames(data.grape)
	CondRank1 = CondRank[, order(CondRank["Med.cond.score", ], decreasing = T)]
	CondRank2 = CondRank[, order(CondRank["Med.cond.score", ], decreasing = F)]
	if(heatmap) {
		Data = as.data.frame(data.grape.scale[j, c(CondRank1[1, 1:Top], CondRank2[1, 1:Top])])
		Data$class = Module
		HeatMap(Data, clusterRow, clusterCol, cutrow, cutcol, address = address, binary = binary)
	}
	A = data.grape.scale[j, c(CondRank1[1, 1:Top])] > 0
	B = data.grape.scale[j, c(CondRank2[1, 1:Top])] > 0
	Score = (sum(A)-sum(B))/length(A)
	RR = data.frame("Score" = Score, t(data.frame(colnames(CondRank1)[1:Top])), t(data.frame(colnames(CondRank2)[1:Top])))
	return(RR)
}


##########################################################################




##########################################################################
# Grape data, expression and the graph
##########################################################################
data.grape = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/Vitis_vinifera/GE.mean.SD.txt", nrows = 23629, comment.char = "", header = T)
cond.devel = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/Vitis_vinifera/development.annotaion.txt")
colnames(data.grape) = cond.devel[,2]
data.grape.scale = Scaling(data.grape)
gene.names = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/Vitis_vinifera/gene.names.txt")
load("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Graph_V2_90.RData")
g.grape = g
##########################################################################



########################################################################## 
# Condition scores
########################################################################## 

list = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Scripts/General/order.condition_score.Grape.txt")
list = as.vector(unlist(list))
## Read text or RData file with read.table or load function
setwd("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Raw/Condition_score/")
# C = lapply(list, read.table, sep= "\t" , header=T)
C = lapply(list, read.table, header = T)
total_cond_scores = matrix(NA, length(unlist(C[1])), length(list))
for (i in 1:length(list)){
	print(i)
	total_cond_scores[, i] = as.matrix(unlist(C[i]))
}
colnames(total_cond_scores) = gene.names[, 1]

########################################################################## 


OrthZhen = read.delim("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Orthology/Grape-Arabidopsis/PollenRootLeafclean.txt")
Score = read.delim("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Compendium_wise/development/PollenRoot90.txt")
Score95 = read.delim("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Compendium_wise/development/PollenRoot95.txt")
PollenLeaf = intersect(Score95[Score95$PollenLeafFDR < 0.05 & Score95$RootLeaf < 10, "Grape"], Score[Score$PollenLeafFDR < 0.05 & Score$RootLeaf < 10, "Grape"])
RootLeaf = intersect(Score95[Score95$RootLeafFDR < 0.01 & Score95$PollenLeaf < 10, "Grape"], Score[Score$RootLeafFDR < 0.01 & Score$PollenLeaf < 10, "Grape"])
candidate = Score[Score$Grape %in% PollenLeaf,]
candidate = Score[Score$Grape %in% RootLeaf,]

candidate["ImageScore"] = NA
SS = foreach(i = 1:nrow(candidate), .combine = rbind) %do% {
	print(i)
	j = data.frame("name" = gene.names[neighborhood(g.grape, 1, candidate$Grape[i], mode = c("out"))[[1]],])
	j = as.character(j[,1])
	jj = as.character(OrthZhen[OrthZhen$Module == "Pollen", "name"])
	jj = as.character(OrthZhen[OrthZhen$Module == "Root", "name"])
	jjj = unique(jj[jj %in% j])
	OrderCondition(jjj, heatmap = T, Module = "Pollen", Top = 15, clusterRow = T, clusterCol = F, cutrow = 1, cutcol = 1, binary = F, address = NA)
}


cond.devel[cond.devel$V1 %in% unique(c(as.matrix(SS[,2:16]))),]


cond.devel[cond.devel$V1 %in% c(as.matrix(SS[1,2:11])),]

candidate = candidate[order(-candidate$ImageScore),]

cond.devel[cond.devel[,1] %in% names(CondRank1[1, 1:10]),]
cond.devel[cond.devel[,1] %in% names(CondRank2[1, 1:10]),]


#####################################################################################################################################################################################
# Heatmap for grapevines orthologues across all conditions
#####################################################################################################################################################################################
library("isa2")
library("foreach")
library("doMC")
library(igraph)
options(width = 200)
library(pheatmap)
library(preprocessCore)

##########################################################################
# Reading Grape expression data and doing the quantile normalization
##########################################################################
data.grape = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/Vitis_vinifera/GE.mean.SD.txt", nrows = 23629, comment.char = "", header = T)
cond.devel = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/Vitis_vinifera/development.annotaion.txt")
colnames(data.grape) = cond.devel[,2]
data.grape.norm = normalize.quantiles(as.matrix(data.grape))
rownames(data.grape.norm) = rownames(data.grape)
colnames(data.grape.norm) = colnames(data.grape)
data.grape.scale = Scaling(data.grape)
data.grape.norm.scale = Scaling(data.grape.norm)
data.grapeEEi = EEi(data.grape)
data.grape.binary = ifelse(data.grape.scale > 0, 1, 0)
data.grape.PW = Scaling(data.grape[, -c(46:54, 85:93, 139:147)])
##########################################################################


##########################################################################
# Remove Berry related conditions
##########################################################################
data.grape.Berry = data.grape[, -c(34:54, 127:147)]
data.grape.Berry.scale = Scaling(data.grape.Berry)
data.grape.Berry.norm = data.grape.norm[, -c(34:54, 127:147)]
data.grape.Berry.scale.norm = Scaling(data.grape.Berry.norm)
data.grape.BerryEEi = EEi(data.grape.Berry)
##########################################################################


data.grape.PollenRoot = Scaling(data.grape[, c(1:3,61:78, 94:96)])
data.grape.PollenRoot.scale = data.grape.scale[, c(1:3,61:78, 94:96)]



##########################################################################
# Removing ORTH with multiple gene names to one PLAZA ID
##########################################################################
OrthZhen = read.delim("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Orthology/Grape-Arabidopsis/PollenRootLeafclean.txt")
PLAZA.Dup = as.character(OrthZhen$PLAZA.ID[duplicated(OrthZhen$PLAZA.ID)])
PLAZA.Dup = PLAZA.Dup[sapply(1:length(PLAZA.Dup), function(i){length(unique(OrthZhen[OrthZhen$PLAZA.ID %in% PLAZA.Dup[i], "name"]))}) > 1]
OrthZhen = OrthZhen[!(OrthZhen$PLAZA.ID %in% PLAZA.Dup),]
map = read.delim("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Orthology/Grape-Arabidopsis/Map.12x.plaza.txt")
map = map[!(duplicated(map[, c("name", "PLAZA.ID")])),]
##########################################################################


##########################################################################
# Extracting APs with corresponding ORTH
##########################################################################
Pollen = unique(as.character(OrthZhen[OrthZhen$Module == "Pollen", "name"]))
Root = unique(as.character(OrthZhen[OrthZhen$Module == "Root", "name"]))
RootPollen = intersect(Pollen, Root)
Pollen = Pollen[!(Pollen %in% RootPollen)]
Root = Root[!(Root %in% RootPollen)]
##########################################################################


##########################################################################
# Marcopolo data
##########################################################################
marcoP = read.delim("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/Marcopolo/Pollen.txt", header = F)
marcoR = read.delim("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/Marcopolo/Root.txt", header = F)
marcoA = read.delim("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/Marcopolo/Root.txt", header = F)
marcoP = read.delim("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/Marcopolo/RootPollen-FlowerRoot.txt", header = F)


MarcoPollen = readMarco(marcoP, "Pollen")
MarcoRoot = readMarco(marcoR, "Root")
MarcoAll = readMarco(marcoA, "All")
Marco = rbind(MarcoPollen, MarcoRoot)

MarcoPollenRoot = readMarco(marcoP, "conserve")
##########################################################################



##########################################################################
# Ancestral expression values, consists of ICC and ranked orthologues
##########################################################################
PolleN = read.delim("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Grape/Multiple_othologous_modules/Pollen.txt")
RooT = read.delim("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Grape/Multiple_othologous_modules/Root.txt")
##########################################################################



ClassPollen = sapply(1:length(Pollen), function(i){as.character(PolleN[PolleN$Grape == Pollen[i], "Module"][1])})
ClassRoot = sapply(1:length(Root), function(i){as.character(RooT[RooT$Grape == Root[i], "Module"][1])})

Data = data.frame(rbind(data.grape.scale[Pollen,], data.grape.scale[Root,]), "class" = c(ClassPollen, ClassRoot))
Data.EEi = data.frame(rbind(data.grapeEEi[Pollen,], data.grapeEEi[Root,]),"class" = c(ClassPollen, ClassRoot))
Data.norm = data.frame(rbind(data.grape.norm.scale[Pollen,], data.grape.norm.scale[Root,]), "class" = c(ClassPollen, ClassRoot))
Data.Berry = data.frame(rbind(data.grape.Berry.scale[Pollen,], data.grape.Berry.scale[Root,]),"class" = c(ClassPollen, ClassRoot))
Data.Berry.norm = data.frame(rbind(data.grape.Berry.scale.norm[Pollen,], data.grape.Berry.scale.norm[Root,]),"class" = c(ClassPollen, ClassRoot))
Data.BerryEEi = data.frame(rbind(data.grape.BerryEEi[Pollen,], data.grape.BerryEEi[Root,]),"class" = c(ClassPollen, ClassRoot))
Data.RF = data.frame(rbind(data.grape.RF.scale[Pollen,], data.grape.RF.scale[Root,]),"class" = c(ClassPollen, ClassRoot))
Data.binary = data.frame(rbind(data.grape.binary[Pollen,], data.grape.binary[Root,]), "class" = c(ClassPollen, ClassRoot))
Data.PollenRoot = data.frame(rbind(data.grape.PollenRoot[Pollen,], data.grape.PollenRoot[Root,]), "class" = c(ClassPollen, ClassRoot))
Data.PollenRoot.scale = data.frame(rbind(data.grape.PollenRoot.scale[Pollen,], data.grape.PollenRoot.scale[Root,]), "class" = c(ClassPollen, ClassRoot))
Data.PW = data.frame(rbind(data.grape.PW[Pollen,], data.grape.PW[Root,]), "class" = c(ClassPollen, ClassRoot))

ClassPollen = sapply(1:length(Pollen), function(i){as.character(PolleN[PolleN$Grape == Pollen[i], "Module"][1])})
DataPollen = data.frame(data.grape.scale[Pollen,], "class" = ClassPollen)
DataPollenEEi = data.frame(data.grapeEEi[Pollen,], "class" = ClassPollen)
DataPollen.Berry = data.frame(data.grape.Berry.scale[Pollen,], "class" = ClassPollen)
DataPollen.BerryEEi = data.frame(data.grape.BerryEEi[Pollen,], "class" = ClassPollen)
DataPollen.Berry.norm = data.frame(data.grape.Berry.scale.norm[Pollen,], "class" = ClassPollen)
DataPollen.PollenRoot = data.frame(data.grape.PollenRoot[Pollen,], "class" = ClassPollen)
DataPollen.PollenRoot.scale = data.frame(data.grape.PollenRoot.scale[Pollen,], "class" = ClassPollen)
DataPollen.RF = data.frame(data.grape.RF.scale[Pollen,], "class" = ClassPollen)
DataPollen.PW = data.frame(data.grape.PW[Pollen,], "class" = ClassPollen)

ClassRoot = sapply(1:length(Root), function(i){as.character(RooT[RooT$Grape == Root[i], "Module"][1])})
DataRoot = data.frame(data.grape.scale[Root,], "class" = ClassRoot)
DataRootEEi = data.frame(data.grapeEEi[Root,], "class" = ClassRoot)
DataRoot.Berry = data.frame(data.grape.Berry.scale[Root,], "class" = ClassRoot)
DataRoot.BerryEEi = data.frame(data.grape.BerryEEi[Root,], "class" = ClassRoot)
DataRoot.Berry.norm = data.frame(data.grape.Berry.scale.norm[Root,], "class" = ClassRoot)
DataRoot.PollenRoot = data.frame(data.grape.PollenRoot[Root,], "class" = ClassRoot)
DataRoot.PollenRoot.scale = data.frame(data.grape.PollenRoot.scale[Root,], "class" = ClassRoot)
DataRoot.RF = data.frame(data.grape.RF.scale[Root,], "class" = ClassRoot)
DataRoot.PW = data.frame(data.grape.PW[Root,], "class" = ClassRoot)


HeatMap(Data, T, T, 3, 5, binary = F, "/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Grape/Plots/Heatmaps/RootPollen.pdf")
HeatMap(Data.EEi, T, T, 3, 5, binary = F,"/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Grape/Plots/Heatmaps/RootPolleEEI.pdf")
HeatMap(Data.norm, T, T, 3, 5, binary = F,"/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Grape/Plots/Heatmaps/RootPolleNorm.pdf")
HeatMap(Data.Berry, T, T, 3, 5, binary = F,"/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Grape/Plots/Heatmaps/RootPolleBerry.pdf")
HeatMap(Data.Berry.norm, T, T, 3, 5, binary = F,"/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Grape/Plots/Heatmaps/RootPolleBerryNorm.pdf")
HeatMap(Data.BerryEEi, T, T, 3, 5, binary = F,"/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Grape/Plots/Heatmaps/RootPolleBerryEEI.pdf")
HeatMap(Data.binary, T, T, 3, 5, binary = T, "/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Grape/Plots/Heatmaps/RootPolleBinary.pdf")
HeatMap(Data.PW, T, T, 3, 5,anntcol = F, binary = F, "/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Grape/Plots/Heatmaps/RootPollenPW.pdf")
HeatMap(Data.PollenRoot, T, T, 3, 5,anntcol = F, binary = F, "/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Grape/Plots/Heatmaps/RootPollen-PollenRoot.pdf")
HeatMap(Data.PollenRoot.scale, T, T, 3, 5, anntcol = F, binary = F, "/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Grape/Plots/Heatmaps/RootPollen-PollenRootScale.pdf")
HeatMap(Data.RF, T, T, 3, 5, binary = F,"/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Grape/Plots/Heatmaps/RootPolleRF.pdf")
HeatMap(MarcoPollenRoot, T, T, 3, 5, anntcol = F, binary = F, "/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Grape/Plots/Heatmaps/RootPollenMarco.pdf")

HeatMap(DataPollen, T, T, 1, 5, anntcol = F,  binary = F,"/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Grape/Plots/Heatmaps/Pollen.pdf")
HeatMap(DataPollenEEi, T, T, 1, 5, binary = F,"/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Grape/Plots/Heatmaps/PollenEEi.pdf")
HeatMap(DataPollen.Berry, T, T, 1, 5, binary = F,"/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Grape/Plots/Heatmaps/PollenBerry.pdf")
HeatMap(DataPollen.BerryEEi, T, T, 1, 5, binary = F,"/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Grape/Plots/Heatmaps/PollenBerryEEi.pdf")
HeatMap(DataPollen.Berry.norm, T, T, 1, 5, binary = F,"/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Grape/Plots/Heatmaps/PollenBerryNorm.pdf")
HeatMap(DataPollen.PW, T, T, 1, 5,anntcol = F, binary = F, "/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Grape/Plots/Heatmaps/PollenPW.pdf")
HeatMap(DataPollen.PollenRoot, T, T, 1, 5,anntcol = F, binary = F, "/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Grape/Plots/Heatmaps/Pollen.PollenRoot.pdf")
HeatMap(DataPollen.PollenRoot.scale, T, T, 1, 5,anntcol = F, binary = F, "/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Grape/Plots/Heatmaps/Pollen.PollenRootScale.pdf")
HeatMap(DataPollen.RF, T, T, 1, 1, binary = F,"/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Grape/Plots/Heatmaps/PollenRF.pdf")
HeatMap(MarcoPollen, T, T, 1, 5, binary = F,"/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Grape/Plots/Heatmaps/PollenMarco.pdf")



HeatMap(DataRoot, T, T, 1, 5, binary = F,"/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Grape/Plots/Heatmaps/Root.pdf")
HeatMap(DataRootEEi, T, T, 1, 5, binary = F,"/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Grape/Plots/Heatmaps/RootEEi.pdf")
HeatMap(DataRoot.Berry, T, T, 1, 5, binary = F,"/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Grape/Plots/Heatmaps/RootBerry.pdf")
HeatMap(DataRoot.BerryEEi, T, T, 1, 5, binary = F,"/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Grape/Plots/Heatmaps/RootBerryEEi.pdf")
HeatMap(DataRoot.Berry.norm, T, T, 1, 5, binary = F,"/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Grape/Plots/Heatmaps/RootBerryNorm.pdf")
HeatMap(DataRoot.PW, T, T, 1, 5,anntcol = F, binary = F, "/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Grape/Plots/Heatmaps/RootPW.pdf")
HeatMap(DataRoot.PollenRoot, T, T, 1, 5,anntcol = F, binary = F, "/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Grape/Plots/Heatmaps/Root.PollenRoot.pdf")
HeatMap(DataRoot.PollenRoot.scale, T, T, 1, 5,anntcol = F, binary = F, "/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Grape/Plots/Heatmaps/Root.PollenRootScale.pdf")
HeatMap(DataRoot.RF, T, T, 1, 1, binary = F,"/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Grape/Plots/Heatmaps/RootRF.pdf")
HeatMap(MarcoRoot, T, T, 1, 5, binary = F,"/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Grape/Plots/Heatmaps/RootMarco.pdf")

##########################################################################





##########################################################################
# Correlation plot for within and between Pollen and Root orthologues
##########################################################################
Data.Corr = data.frame(cor(t(data.grape.scale[c(Pollen, Root),])), "class" = c(ClassPollen, ClassRoot))
Data.CorrPollen = data.frame(cor(t(data.grape.scale[c(Pollen),])), "class" = ClassPollen)
Data.CorrRoot = data.frame(cor(t(data.grape.scale[c(Root),])), "class" = ClassRoot)

HeatMap(Data.Corr, T, T, 4, 4, binary = F, anntcol=T, address = "/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Grape/Plots/Heatmaps/RootPollenCorrelation.pdf")
HeatMap(Data.CorrPollen, T, T, 1, 1, binary = F, anntcol=T, address = "/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Grape/Plots/Heatmaps/PollenCorrelation.pdf")
HeatMap(Data.CorrRoot, T, T, 1, 1, binary = F, anntcol=T, address = "/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Grape/Plots/Heatmaps/RootCorrelation.pdf")

#####################################################################################################################################################################################
