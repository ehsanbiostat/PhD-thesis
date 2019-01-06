

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

source("/ngsprojects/iwt/ehsab/Gene_duplication/Scripts/General/myImagePlot_sym.R")
source("/ngsprojects/iwt/ehsab/Gene_duplication/Scripts/General/iwanthue.R")

Similarcondition = function(j, CondRank0) {
	TopLeaf = foreach(k = 1:10, .combine = rbind) %do% {
		print(k) 
		AA = foreach(i = 1:ncol(data.exp.scale.dev), .combine = rbind) %do% { 
			data.frame("Cor" = cor(data.exp.scale[j,CondRank0[1,k]], data.exp.scale.dev[j, i]), "Dist" = dist(t(cbind(data.exp.scale[j,CondRank0[1,k]], data.exp.scale.dev[j, i])))[1],
			 "QueryCondition" = colnames(CondRank0)[k], "Condition" = colnames(data.exp.scale.dev)[i])
		}
		AA
	}
	GSEA = cbind(aggregate(Cor ~ Condition, data = TopLeaf, mean), aggregate(Cor ~ Condition, data = TopLeaf, sd)[,2])
	GSEA = GSEA[order(-GSEA[,2]),]
	return(GSEA)
}


Annot = function(x) {
	Annot = foreach(i = 1:length(x), .combine = c) %do% {
		strsplit(substring(x, 2), "[.]")[[i]][1]
	}
	return(Annot)
}



ReorderMatrix = function(mydata, condReverse, Row, Column) { 
	# set the custom distance and clustering functions, per your example
	hclustfunc <- function(x) hclust(x, method="ward.D2")
	distfunc <- function(x) dist(x, method="euclidean")

	# perform clustering on rows and columns
	cl.row <- hclustfunc(distfunc(t(mydata)))
	cl.col <- hclustfunc(distfunc(mydata))
	if(Row & Column & condReverse) mydata = mydata[cl.col$order, rev(cl.row$order)]
	if(!Row & Column & !condReverse) mydata = mydata[, cl.row$order]
	if(!Row & Column & condReverse) mydata = mydata[, rev(cl.row$order)]
	if(Row & !Column & !condReverse) mydata = mydata[cl.col$order, ]
	if(Row & Column & !condReverse) mydata = mydata[cl.col$order, cl.row$order]
	if(Row & !Column & condReverse) mydata = mydata[rev(cl.col$order), ]
	return(mydata)
}



TopConditions = function(memb.comun, memb.comun.para, TOP) {
	CondRank1 = foreach(i = 1:length(memb.comun), .combine = rbind) %do%{
		rank(total_cond_scores[, memb.comun[i]])
	}
	CondRank2 = foreach(i = 1:length(memb.comun.para), .combine = rbind) %do%{
		rank(total_cond_scores[, memb.comun.para[i]])
	}
			
	
	Med.cond.score1 = apply(CondRank1, 2, median)
	Med.cond.score2 = apply(CondRank2, 2, median)

	memb.comun.para.cond = colnames(data.exp.scale)[which(rank(Med.cond.score1 - Med.cond.score2, ties.method= "random") < TOP+1)]
	memb.comun.cond = colnames(data.exp.scale)[which(rank(Med.cond.score1 - Med.cond.score2, ties.method= "random") > ncol(data.exp.scale)-TOP)]
	Final = data.frame(memb.comun.para.cond, memb.comun.cond)
	return(Final)
}


RandomForest = function(data, x, y, namex, namey, Top, address, classes, title) {
	if(classes == 2) Data = data.frame(rbind.data.frame(data[x,], data[y,]), "class" = c(rep(namex, length(x)), rep(namey, length(y))))
	if(classes > 2) Data = data
	Result = randomForest(factor(class) ~ ., data = Data, ntrees = 20000)
	if(classes == 2) Var = data.frame(colnames(data), varUsed(Result))
	if(classes > 2) Var = data.frame(colnames(data[, -ncol(data)]), varUsed(Result))
	Var = Var[order(-Var[,2]),]
	# Erro = PollenLeaf$err.rate[500,1]
	# myImagePlot_sym(Data[,as.numeric(rownames(Var[1:Top,]))])
	annotationRow = data.frame("Module" = Data$class)
	ColorRamp = rgb(seq(0,1,length=256),seq(0,1,length=256),seq(1,0,length=256))
	ColorLevels = seq(-max(abs(Data[,-ncol(Data)])), max(abs(Data[,-ncol(Data)])), length=length(ColorRamp)+1)
	rownames(annotationRow) = rownames(Data)
	selected = Data[,as.numeric(rownames(Var[1:Top,]))]
	Col.n = GSM = c()
	for(i in 1:Top) {
		col.n = COR.annot.PO[names(COR.annot.PO) %in% Annot(colnames(selected)[i])][[1]]
		Col.n[i] = ifelse(length(col.n) == 1, col.n,
			ifelse(length(col.n) == 2, paste(col.n[1], col.n[2], sep = " | "),
				ifelse(length(col.n) == 3, paste(col.n[1], col.n[2], col.n[3], sep = " | "), colnames(selected)[i])))
		GSM[i] = COR.annot.GSM[names(COR.annot.GSM) %in% Annot(colnames(selected)[i])][[1]]
	}
	# RR = data.frame(Col.n, "A" = unlist(COR.annot.GSM[names(COR.annot.GSM) %in% Annot(colnames(selected))]))
	RR = data.frame(Col.n, GSM)
	RR = paste(RR[,1], RR[,2], sep = " * ")
	colnames(selected) = RR
	pheatmap(selected, cluster_row = F, clustering_method="ward.D2", annotation_row = annotationRow, border_color=NA, cluster_col = T,
		color = ColorRamp, breaks = ColorLevels, treeheight_col = 0, filename = address, fontsize = 5, main = title,
		cellwidth = 3, cellheight = 5, fontsize_col = 2, fontsize_row = 5,  gaps_row = cumsum(table(factor(annotationRow[,1], levels = as.character(unique(annotationRow[,1]))))))
	return(selected)
}


cellheight = 0.17



RandomForestBoost = function(data, x, y, namex, namey, ntrees, classes) {
	require(varSelRF)
	require(parallel)
	forkCL <- makeForkCluster(4)
	clusterSetRNGStream(forkCL, iseed = 123)
	clusterEvalQ(forkCL, library(varSelRF))
	if(classes == 2) Data = rbind.data.frame(data[x,], data[y,])
	if(classes > 2) Data = data
	Class = as.factor(c(rep(namex, length(x)), rep(namey, length(y))))
	Result = varSelRFBoot(Data, Class, ntrees = ntrees, TheCluster = forkCL)
	return(Result)
	stopCluster(forkCL)
}


require(parallel)
forkCL <- makeForkCluster(2)
clusterSetRNGStream(forkCL, iseed = 123)
clusterEvalQ(forkCL, library(varSelRF))
	
RandomForestSel = function(data, x, y, namex, namey, ntrees, classes, sd) {
	require(varSelRF)
	if(classes == 2) Data = rbind.data.frame(data[x,], data[y,])
	if(classes > 2) Data = data
	Class = as.factor(c(rep(namex, length(x)), rep(namey, length(y))))
	Result = varSelRF(Data, Class, ntree = ntrees, c.sd = sd, whole.range = FALSE)
	return(Result)
}

stopCluster(forkCL)



CirculHeatmap = function(GenExp, Reorder, condReverse, Module, address) {
	mat = as.matrix(GenExp)
	if(Reorder) mat = as.matrix(ReorderMatrix(GenExp, condReverse, Row = F, Column = T)) # clustered matrixed on and column
	mat = mat[c(seq(1,nrow(mat)/2),seq(nrow(mat),(nrow(mat)/2)+1)),]
	ColorRamp = rgb(seq(0,1,length=256),seq(0,1,length=256),seq(1,0,length=256))
	ColorLevels = seq(-max(mat), max(mat), length=length(ColorRamp))
	col_fun = colorRamp2(ColorLevels, ColorRamp)
	col_mat = col_fun(mat)
	color = iwanthue(nrow(mat)/2)
	circos.clear()
	pdf(address, width = 10, height = 10)
	par(mar = c(1, 1, 1, 1))
	circos.par(cell.padding = c(0, 0, 0, 0), gap.degree = 1, start.degree = 90, clock.wise = TRUE)
	factors = letters[1]
	circos.par(points.overflow.warning = FALSE)
	circos.initialize(factors = factors, xlim = c(0, nrow(mat)))
	# Drawing the heatmap
	circos.trackPlotRegion(factors = NULL, ylim = c(0, ncol(mat)), track.height = 0.7, bg.border = NA, panel.fun = function(x, y) {
		nr = nrow(mat)
		nc = ncol(mat)
		for(i in 1:nc) {
	    	for(j in 1:nr) {
	        	circos.rect(j-1, nc-i, j, nc-i+1, border = col_mat[j, i], col = col_mat[j, i])
	    	}
		}
	})
	if(Module == "Pollen") {
		for (i in 1:c(nrow(mat)/2)) { # Anchorpoint connections in the middle of plot
			circos.link("a", which(rownames(mat) %in% as.character(AnchPollen$Pollen[i])) - 0.5,
				"a", which(rownames(mat) %in% as.character(AnchPollen$Leaf[i])) - 0.5, rou = 0.29, col = color[i])
		}
	}
	if(Module == "Root") {
		for (i in 1:c(nrow(mat)/2)) { # Anchorpoint connections in the middle of plot
			circos.link("a", which(rownames(mat) %in% as.character(AnchRoot$Root[i])) - 0.5,
				"a", which(rownames(mat) %in% as.character(AnchRoot$Leaf[i])) - 0.5, rou = 0.29, col = color[i])
		}
	}
	circos.axis(# Writting gene names on the outside of the plot
		sector.index = "a", labels = rownames(mat),	major.at = seq(0.5, nrow(mat), by = 1),
		labels.cex = 0.25, major.tick = FALSE,  direction = "outside", labels.facing = "outside",
		labels.away.percentage = 0, major.tick.percentage = 0.01
	)
	ylabels = sapply(1:length(rev(colnames(mat))), function(i) {strsplit(rev(colnames(mat))[i], "[*]")[[1]][1]})
	for(i in 1:length(ylabels)) { #y axis labels, in this case it is the condition names
		circos.text(0,i-0.5, ylabels[i], sector = "a", facing = "bending.inside", cex = 0.4, col = "black")
	}
	dev.off()
}





#################################################################################################

k = 6
condition.names = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/condition_name.txt", header = T)
NAME = names(table(condition.names))
NAME.spe = NAME[k]
source("/ngsprojects/iwt/ehsab/Gene_duplication/Scripts/General/myImagePlot.R")
source("/ngsprojects/iwt/ehsab/Gene_duplication/Scripts/General/myImagePlot_sym.R")
ref = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/Reference_SD.txt", header = T)

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


# Sig.all = ifelse(data.exp.scale.all > 0, 1, ifelse(data.exp.scale.all < 0, -1, 0))
CORNET.annot = read.delim("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/cornet_allMA_desc_120509.txt")
load("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/PO_CORNET_Annotation.RData")
load("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/GSM_CORNET_Annotation.RData")

# Anchorpoints within two modules
AnchPollen = read.delim("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/Anchorpoint/Geneset1_2_53_90.txt.txt")
AnchRoot = read.delim("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/Anchorpoint/AT1G50090_AT3G19710_52_90.txt")
AnchRootPollen = read.delim("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/Anchorpoint/AnchRootPollen_90.txt")
AnchPollenRef = AnchPollen$Ref
AnchPollenLeafRef = AnchPollen$References
AnchRootRef = AnchRoot$Ref
AnchRootLeafRef = AnchRoot$References


AnchRootPollen = intersect(AnchRoot$Root, AnchPollen$Pollen)
 AnchLeafLeaf = intersect(AnchRoot$Leaf, AnchPollen$Leaf)
 AnchRootS = unique(as.character(AnchRoot$Root[!(AnchRoot$Root %in% AnchRootPollen)]))
 AnchPollenS = unique(as.character(AnchPollen$Pollen[!(AnchPollen$Pollen %in% AnchRootPollen)]))
 AnchRootLeafS = unique(as.character(AnchRoot$Leaf[!(AnchRoot$Leaf %in% AnchLeafLeaf)]))
 AnchPollenLeafS = unique(as.character(AnchPollen$Leaf[!(AnchPollen$Leaf %in% AnchLeafLeaf)]))



# Top conditions based on the most distance condition scores for APs
PollenLeafTop10 = TopConditions(AnchPollenRef, AnchPollenLeafRef) 
RootLeafTop10 = TopConditions(AnchRootRef, AnchRootLeafRef) 

write.table(CORNET.annot[CORNET.annot$Name %in% Annot(as.vector(RootLeafTop10[,1])),], file = "/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Top10_RootLeaf.txt", row.names = F, quote = F, sep = "\t")
myImagePlot_sym(data.exp.scale[c(AnchRootRef, AnchRootLeafRef), c(which(colnames(data.exp.scale) %in% colnames(data.exp.scale)[which(rank(Med.cond.score1, ties.method="random") > 124)]), which(colnames(data.exp.scale) %in% RootLeafTop10[,1]))])



colnames(AnchPollen)[2] = "Root/Pollen"
colnames(AnchRoot)[4] = "Root/Pollen"
Anch = rbind(AnchPollen, AnchRoot)
Anch = Anch[!(duplicated(Anch)), ]

AnchTop10 = TopConditions(Anch$Ref, Anch$References) 
myImagePlot_sym(data.exp.scale[c(Anch$Ref, Anch$References), c(which(colnames(data.exp.scale) %in% AnchTop10[,2]), which(colnames(data.exp.scale) %in% AnchTop10[,1]))])



#######################################################################################################################################################################
# Classification method to select Top conditions
#######################################################################################################################################################################


PollenLeaf = RandomForest(classes = 2, data.exp.scale, as.character(AnchPollen$Pollen), as.character(AnchPollen$Leaf), "Pollen", "Leaf", 50,
	title = "Pollen - Leaf modules", 
	address = "/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Plot/Anchorpoint/RFTop10/TopConditionsPollenRF1.pdf")
PollenLeafBoost = RandomForestBoost(classes = 2, data.exp.scale, as.character(AnchPollen$Pollen), as.character(AnchPollen$Leaf), "Pollen", "Leaf")

RootLeaf = RandomForest(classes = 2, data.exp.scale, as.character(AnchRoot$Root), as.character(AnchRoot$Leaf), "Root", "Leaf", 50,
	title = "Root - Leaf modules", 
	address = "/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Plot/Anchorpoint/RFTop10/TopConditionsRootRF1.pdf")
RootLeafBoost = RandomForestBoost(classes = 2, data.exp.scale, as.character(AnchRoot$Root), as.character(AnchRoot$Leaf), "Root", "Leaf")


RootPollen = RandomForest(classes = 2, data.exp.scale, as.character(AnchRootPollen$Root.Pollen), as.character(AnchRootPollen$Leaf), "Root & Pollen", "Leaf", 50,
	title = "RootPollen - Leaf modules", 
	address = "/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Plot/Anchorpoint/RFTop10/TopConditionsRootPollenRF1.pdf")


LeafLeaf = RandomForest(classes = 2, data.exp.scale, as.character(AnchPollen$Leaf), as.character(AnchRoot$Leaf), "PollenLeaf", "RootLeaf", 10,
	title = "Leaf - Leaf modules", 
	address = "/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Plot/Anchorpoint/RFTop10/TopConditionsLeafLeafRF1.pdf")

LeafLeafBoost = RandomForestBoost(classes = 2, data.exp.scale, as.character(AnchPollen$Leaf), as.character(AnchRoot$Leaf), "PollenLeaf", "RootLeaf")


PollenLeafAdd = RandomForest(classes = 2, data.exp.scale.dev, as.character(AnchPollen$Pollen), as.character(AnchPollen$Leaf), "Pollen", "Leaf", 30,
	title = "Divergence pattern for the Pollen and Leaf modules",
	address = "/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Plot/Anchorpoint/RFaddCond/SimilarConditionsPollenRF.pdf")
RootLeafAdd = RandomForest(classes = 2, data.exp.scale.dev, as.character(AnchRoot$Root), as.character(AnchRoot$Leaf), "Root", "Leaf", 30,
	title = "Divergence pattern for the Root and Leaf modules",
	address = "/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Plot/Anchorpoint/RFaddCond/SimilarConditionsRootRF.pdf")
RootPollenAdd = RandomForest(classes = 2, data.exp.scale.dev, as.character(AnchRootPollen$Root.Pollen), as.character(AnchRootPollen$Leaf), "Root & Pollen", "Leaf", 30,
	title = "Divergence pattern for the Root & Pollen and Leaf modules",
	address = "/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Plot/Anchorpoint/RFaddCond/SimilarConditionsRootPollenRF.pdf")

# Gene coexpression divergence pattern for 

Data = data.frame(rbind(data.exp.scale[AnchRootPollen,],
	data.exp.scale[AnchPollenS,],data.exp.scale[AnchRootS,], data.exp.scale[AnchLeafLeaf,],
	data.exp.scale[AnchPollenLeafS,],data.exp.scale[AnchRootLeafS,]),
"class" = c(rep("RootPollen", length(AnchRootPollen)),
		rep("Pollen", length(AnchPollenS)),
		rep("Root", length(AnchRootS)),
		rep("LeafLeaf", length(AnchLeafLeaf)),
		rep("PollenLeaf", length(AnchPollenLeafS)),
		rep("RootLeaf", length(AnchRootLeafS))))

PollenRoot = RandomForest(classes = 4, data = Data, title = "", Top = 40, address = "/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Plot/Anchorpoint/RFTop10/TopConditionsPollenRootRF.pdf")

Data.add = data.frame(rbind(data.exp.scale.dev[AnchRootPollen,],
	data.exp.scale.dev[AnchPollenS,],data.exp.scale.dev[AnchRootS,], data.exp.scale.dev[AnchLeafLeaf,],
	data.exp.scale.dev[AnchPollenLeafS,],data.exp.scale.dev[AnchRootLeafS,]),
"class" = c(rep("RootPollen", length(AnchRootPollen)),
		rep("Pollen", length(AnchPollenS)),
		rep("Root", length(AnchRootS)),
		rep("LeafLeaf", length(AnchLeafLeaf)),
		rep("PollenLeaf", length(AnchPollenLeafS)),
		rep("RootLeaf", length(AnchRootLeafS))))



Data.add = data.frame(rbind(data.exp.scale.dev[AnchRootPollen,], data.exp.scale.dev[AnchLeafLeaf,],
	data.exp.scale.dev[AnchPollenLeafS,],data.exp.scale.dev[AnchRootS,],
	data.exp.scale.dev[AnchPollenS,],data.exp.scale.dev[AnchRootLeafS,]),
"class" = c(rep("RootPollen", length(AnchRootPollen)),
		rep("LeafLeaf", length(AnchLeafLeaf)),
		rep("PollenLeaf", length(AnchPollenLeafS)),
		rep("Root", length(AnchRootS)),
		rep("Pollen", length(AnchPollenS)),
		rep("RootLeaf", length(AnchRootLeafS))))

PollenRootAdd = RandomForest(classes = 4, data = Data.add, Top = 40, address = "/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Plot/Anchorpoint/RFaddCond/SimilarConditionsPollenRootRF.pdf")

#######################################################################################################################################################################



#######################################################################################################################################################################
# Circular Heatmap
#
#######################################################################################################################################################################

AnchPollen$Pollen = sapply(1:nrow(AnchPollen), function(i) {ifelse(duplicated(AnchPollen$Pollen)[i], paste(AnchPollen$Pollen[i],1,sep = ""), as.character(AnchPollen$Pollen[i]))})
AnchPollen$Leaf = sapply(1:nrow(AnchPollen), function(i) {ifelse(duplicated(AnchPollen$Leaf)[i], paste(AnchPollen$Leaf[i],1,sep = ""), as.character(AnchPollen$Leaf[i]))})

AnchRoot$Root = sapply(1:nrow(AnchRoot), function(i) {ifelse(duplicated(AnchRoot$Root)[i], paste(AnchRoot$Root[i],1,sep = ""), as.character(AnchRoot$Root[i]))})
AnchRoot$Leaf = sapply(1:nrow(AnchRoot), function(i) {ifelse(duplicated(AnchRoot$Leaf)[i], paste(AnchRoot$Leaf[i],1,sep = ""), as.character(AnchRoot$Leaf[i]))})



CirculHeatmap(PollenLeaf, Reorder = T, condReverse = T, Module = "Pollen",
	"/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Plot/Anchorpoint/RFTop10/TopConditionsPollenRFCirclularClustered.pdf")
CirculHeatmap(PollenLeaf, Reorder = T, condReverse = F, Module = "Pollen",
	"/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Plot/Anchorpoint/RFTop10/TopConditionsPollenRFCirclular.pdf")

CirculHeatmap(RootLeaf, Reorder = T, condReverse = T, Module = "Root",
	"/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Plot/Anchorpoint/RFTop10/TopConditionsRootRFCirclularClustered.pdf")
CirculHeatmap(RootLeaf, Reorder = T, condReverse = F, Module = "Root",
	"/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Plot/Anchorpoint/RFTop10/TopConditionsRootRFCirclular.pdf")



CirculHeatmap(PollenLeafAdd, Reorder = T, condReverse = T, Module = "Pollen",
	"/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Plot/Anchorpoint/RFaddCond/SimilarConditionsPollenRFCirclular.pdf")
CirculHeatmap(RootLeafAdd, Reorder = T, condReverse = F, Module = "Root",
	"/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Plot/Anchorpoint/RFaddCond/SimilarConditionsRootRFCirclular.pdf")




#######################################################################################################################################################################
# Arcplot for heatmaps
#######################################################################################################################################################################

AnchPollen = read.delim("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/Anchorpoint/Geneset1_2_53_90.txt.txt")
AnchRoot = read.delim("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/Anchorpoint/AT1G50090_AT3G19710_52_90.txt")

edgelistPollen = as.matrix(AnchPollen[, c(2,4)])
PollenNodes = unique(c(as.matrix(AnchPollen[, c(2,4)] )))
new_ord = sample(1:length(PollenNodes), length(PollenNodes))

pdf("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Plot/Anchorpoint/Arcs/Pollen.pdf", width = 3, height = 7)
arcplot(edgelistPollen, vertices = PollenNodes, show.labels = FALSE, horizontal= F, show.nodes=F, outer = T)
dev.off()

png("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Plot/Anchorpoint/Arcs/Pollen.png", width = 5, height = 10, units = "in", res = 200)
arcplot(edgelistPollen, vertices = PollenNodes, show.labels = FALSE, horizontal= F, show.nodes=F, outer = T)
dev.off()

edgelistRoot = as.matrix(AnchRoot[, c(2,4)])
RootNodes = unique(c(as.matrix(AnchRoot[, c(2,4)] )))

pdf("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Plot/Anchorpoint/Arcs/Root.pdf", width = 3, height = 7)
arcplot(edgelistRoot, vertices = RootNodes, show.labels = FALSE, horizontal= F, show.nodes=F, outer = T)
dev.off()

png("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Plot/Anchorpoint/Arcs/Root.png", width = 5, height = 10, units = "in", res = 200)
arcplot(edgelistRoot, vertices = RootNodes, show.labels = FALSE, horizontal= F, show.nodes=F, outer = T)
dev.off()



# Rotate and crop the png file
https://pixlr.com/editor/

# Merge heatmap and arcplot
pdfjam \
"Drives/research/deepseq/scratch/scratch-iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Plot/Anchorpoint/Arcs/Root.pdf" \
"Drives/research/deepseq/scratch/scratch-iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Plot/Anchorpoint/RFTop10/TopConditionsRootRF.pdf" \
--nup 2x1 --landscape  --outfile "Drives/research/deepseq/scratch/scratch-iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Plot/Anchorpoint/Arcs/RootArc.pdf"
pdfjam \
"Drives/research/deepseq/scratch/scratch-iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Plot/Anchorpoint/Arcs/Pollen.pdf" \
"Drives/research/deepseq/scratch/scratch-iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Plot/Anchorpoint/RFTop10/TopConditionsPollenRF.pdf" \
--nup 2x1 --landscape  --outfile "Drives/research/deepseq/scratch/scratch-iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Plot/Anchorpoint/Arcs/PollenArc.pdf"


# Move the arcplot close enough to the heatmap
https://www.pdfescape.com/

# Merge Pollen and Root
pdfjam \
"Drives/research/deepseq/scratch/scratch-iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Plot/Anchorpoint/Arcs/PollenArc.pdf" \
"Drives/research/deepseq/scratch/scratch-iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Plot/Anchorpoint/Arcs/RootArc.pdf" \
--nup 2x1 --landscape  --outfile "Drives/research/deepseq/scratch/scratch-iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Plot/Anchorpoint/Arcs/PollenRootArc.pdf"
pdfjam \
"Drives/research/deepseq/scratch/scratch-iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Plot/Anchorpoint/Arcs/PollenArc.pdf" \
"Drives/research/deepseq/scratch/scratch-iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Plot/Anchorpoint/Arcs/RootArc.pdf" \
--nup 2x1 --landscape  --outfile "Drives/research/deepseq/scratch/scratch-iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Plot/Anchorpoint/Arcs/PollenRootArc.pdf"

########################################################################################################################################################################################






#######################################################################################################################################################################################
# For entire module (Pollen-Leaf & Root-Leaf)
#######################################################################################################################################################################################
Pollen = read.table("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/Module/Pollen_90.txt")
PollenLeaf = read.table("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/Module/Leaf_90.txt")
Root = read.table("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/Module/AT1G50090_Root_90.txt")
RootLeaf = read.table("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/Module/AT3G19710_Leaf_90.txt")


pollen = RandomForest(classes = 2, data.exp.scale, as.character(Pollen[,1]), as.character(PollenLeaf[,1]), "Pollen", "Leaf", 30,
title = "Pollen - Root modules", 
address = "/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Plot/Module/Top10/TopConditionsPollenRootRF.pdf")


root = RandomForest(classes = 2, data.exp.scale, as.character(Root[,1]), as.character(RootLeaf[,1]), "Root", "Leaf", 30,
title = "Root - Leaf modules", 
address = "/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Plot/Module/Top10/TopConditionsRootRF.pdf")


leaf = RandomForest(classes = 2, data.exp.scale, as.character(PollenLeaf[,1]), as.character(RootLeaf[,1]), "PollenLeaf", "RootLeaf", 10,
title = "PollenLeafRoot - RootLeaf modules", 
address = "/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Plot/Module/Top10/TopConditionsLeafLeaf.pdf")


RootLeafS = setdiff(as.character(RootLeaf[,1]), as.character(PollenLeaf[,1]))
PollenLeafS = setdiff(as.character(PollenLeaf[,1]), as.character(RootLeaf[,1]))

leafS = RandomForest(classes = 2, data.exp.scale, PollenLeafS, RootLeafS, "PollenLeaf", "RootLeaf", 10,
title = "PollenLeaf - RootLeaf modules Specific", 
address = "/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Plot/Module/Top10/TopConditionsLeafLeafSpecific.pdf")

#######################################################################################################################################################################################





#######################################################################################################################################################################################
# All possible combitation for four modules
#######################################################################################################################################################################################


# overlapped genes
PollenRootOverlap = intersect(Pollen, Root)
PollenRootLeafOverlap = intersect(Pollen, RootLeaf)
RootPollenLeafOverlap = intersect(PollenLeaf, Root)
PollenLeafRootLeafOverlap = intersect(PollenLeaf, RootLeaf)

PollenS = setdiff(Pollen, PollenRootOverlap)
RootS = setdiff(Root, PollenRootOverlap)
PollenLeafS = setdiff(PollenLeaf, PollenLeafRootLeafOverlap)
RootLeafS = setdiff(RootLeaf, PollenLeafRootLeafOverlap)


# Shared duplicated genes
PollenDup = as.character(unique(GeneFamilyDup[GeneFamilyDup$Var1 %in% Pollen,"Var1"]))
PollenLeafDup = as.character(unique(GeneFamilyDup[GeneFamilyDup$Var1 %in% PollenLeaf,"Var1"]))
RootDup = as.character(unique(GeneFamilyDup[GeneFamilyDup$Var1 %in% Root,"Var1"]))
RootLeafDup = as.character(unique(GeneFamilyDup[GeneFamilyDup$Var1 %in% RootLeaf,"Var1"]))
PollenPollen = intersect(as.character(unique(GeneFamilyDup[GeneFamilyDup$Var1 %in% Pollen, "Var2"])), Pollen)
RootRoot = intersect(as.character(unique(GeneFamilyDup[GeneFamilyDup$Var1 %in% Root, "Var2"])), Root)
LeafLeafPollen = intersect(as.character(unique(GeneFamilyDup[GeneFamilyDup$Var1 %in% PollenLeaf, "Var2"])), PollenLeaf)
LeafLeafRoot = intersect(as.character(unique(GeneFamilyDup[GeneFamilyDup$Var1 %in% RootLeaf, "Var2"])), RootLeaf)
RootPollen = intersect(as.character(unique(GeneFamilyDup[GeneFamilyDup$Var1 %in% Pollen, "Var2"])), Root)
PollenRoot = intersect(as.character(unique(GeneFamilyDup[GeneFamilyDup$Var1 %in% Root, "Var2"])), Pollen)
LeafPollen = intersect(as.character(unique(GeneFamilyDup[GeneFamilyDup$Var1 %in% Pollen, "Var2"])), PollenLeaf)
PollenRootLeaf = intersect(as.character(unique(GeneFamilyDup[GeneFamilyDup$Var1 %in% RootLeaf, "Var2"])), Pollen)
RootLeafPollen = intersect(as.character(unique(GeneFamilyDup[GeneFamilyDup$Var1 %in% Pollen, "Var2"])), RootLeaf)
Pollenleaf = intersect(as.character(unique(GeneFamilyDup[GeneFamilyDup$Var1 %in% PollenLeaf, "Var2"])), Pollen)
LeafRoot = intersect(as.character(unique(GeneFamilyDup[GeneFamilyDup$Var1 %in% Root, "Var2"])), RootLeaf)
Rootleaf = intersect(as.character(unique(GeneFamilyDup[GeneFamilyDup$Var1 %in% RootLeaf, "Var2"])), Root)
RootPollenLeaf = intersect(as.character(unique(GeneFamilyDup[GeneFamilyDup$Var1 %in% PollenLeaf, "Var2"])), Root)
PollenLeafRoot = intersect(as.character(unique(GeneFamilyDup[GeneFamilyDup$Var1 %in% Root, "Var2"])), PollenLeaf)
RootLeafPollenLeaf = intersect(as.character(unique(GeneFamilyDup[GeneFamilyDup$Var1 %in% PollenLeaf, "Var2"])), RootLeaf)
PollenLeafRootLeaf = intersect(as.character(unique(GeneFamilyDup[GeneFamilyDup$Var1 %in% RootLeaf, "Var2"])), PollenLeaf)


# Shared duplicated genes Large Scale Duplication (LSD)
PollenDupLSD = as.character(unique(paralogous[paralogous$Ref %in% Pollen,"Ref"]))
PollenLeafDupLSD = as.character(unique(paralogous[paralogous$Ref %in% PollenLeaf,"Ref"]))
RootDupLSD = as.character(unique(paralogous[paralogous$Ref %in% Root,"Ref"]))
RootLeafDupLSD = as.character(unique(paralogous[paralogous$Ref %in% RootLeaf,"Ref"]))
PollenPollenLSD = intersect(as.character(unique(paralogous[paralogous$Ref %in% Pollen, "Anch"])), Pollen)
RootRootLSD = intersect(as.character(unique(paralogous[paralogous$Ref %in% Root, "Anch"])), Root)
LeafLeafPollenLSD = intersect(as.character(unique(paralogous[paralogous$Ref %in% PollenLeaf, "Anch"])), PollenLeaf)
LeafLeafRootLSD = intersect(as.character(unique(paralogous[paralogous$Ref %in% RootLeaf, "Anch"])), RootLeaf)
RootPollenLSD = intersect(as.character(unique(paralogous[paralogous$Ref %in% Pollen, "Anch"])), Root)
PollenRootLSD = intersect(as.character(unique(paralogous[paralogous$Ref %in% Root, "Anch"])), Pollen)
LeafPollenLSD = intersect(as.character(unique(paralogous[paralogous$Ref %in% Pollen, "Anch"])), PollenLeaf)
PollenRootLeafLSD = intersect(as.character(unique(paralogous[paralogous$Ref %in% RootLeaf, "Anch"])), Pollen)
RootLeafPollenLSD = intersect(as.character(unique(paralogous[paralogous$Ref %in% Pollen, "Anch"])), RootLeaf)
PollenleafLSD = intersect(as.character(unique(paralogous[paralogous$Ref %in% PollenLeaf, "Anch"])), Pollen)
LeafRootLSD = intersect(as.character(unique(paralogous[paralogous$Ref %in% Root, "Anch"])), RootLeaf)
RootleafLSD = intersect(as.character(unique(paralogous[paralogous$Ref %in% RootLeaf, "Anch"])), Root)
RootPollenLeafLSD = intersect(as.character(unique(paralogous[paralogous$Ref %in% PollenLeaf, "Anch"])), Root)
PollenLeafRootLSD = intersect(as.character(unique(paralogous[paralogous$Ref %in% Root, "Anch"])), PollenLeaf)
RootLeafPollenLeafLSD = intersect(as.character(unique(paralogous[paralogous$Ref %in% PollenLeaf, "Anch"])), RootLeaf)
PollenLeafRootLeafLSD = intersect(as.character(unique(paralogous[paralogous$Ref %in% RootLeaf, "Anch"])), PollenLeaf)




# Shared duplicated genes with no overlap
PollenDupS = as.character(unique(GeneFamilyDup[GeneFamilyDup$Var1 %in% PollenS,"Var1"]))
PollenLeafDupS = as.character(unique(GeneFamilyDup[GeneFamilyDup$Var1 %in% PollenLeafS,"Var1"]))
RootDupS = as.character(unique(GeneFamilyDup[GeneFamilyDup$Var1 %in% RootS,"Var1"]))
RootLeafDupS = as.character(unique(GeneFamilyDup[GeneFamilyDup$Var1 %in% RootLeafS,"Var1"]))
PollenPollenS = intersect(as.character(unique(GeneFamilyDup[GeneFamilyDup$Var1 %in% PollenS, "Var2"])), PollenS)
RootRootS = intersect(as.character(unique(GeneFamilyDup[GeneFamilyDup$Var1 %in% RootS, "Var2"])), RootS)
LeafLeafPollenS = intersect(as.character(unique(GeneFamilyDup[GeneFamilyDup$Var1 %in% PollenLeafS, "Var2"])), PollenLeafS)
LeafLeafRootS = intersect(as.character(unique(GeneFamilyDup[GeneFamilyDup$Var1 %in% RootLeafS, "Var2"])), RootLeafS)
RootPollenS = intersect(as.character(unique(GeneFamilyDup[GeneFamilyDup$Var1 %in% PollenS, "Var2"])), RootS)
PollenRootS = intersect(as.character(unique(GeneFamilyDup[GeneFamilyDup$Var1 %in% RootS, "Var2"])), PollenS)
LeafPollenS = intersect(as.character(unique(GeneFamilyDup[GeneFamilyDup$Var1 %in% PollenS, "Var2"])), PollenLeafS)
PollenRootLeafS = intersect(as.character(unique(GeneFamilyDup[GeneFamilyDup$Var1 %in% RootLeafS, "Var2"])), PollenS)
RootLeafPollenS = intersect(as.character(unique(GeneFamilyDup[GeneFamilyDup$Var1 %in% PollenS, "Var2"])), RootLeafS)
PollenleafS = intersect(as.character(unique(GeneFamilyDup[GeneFamilyDup$Var1 %in% PollenLeafS, "Var2"])), PollenS)
LeafRootS = intersect(as.character(unique(GeneFamilyDup[GeneFamilyDup$Var1 %in% RootS, "Var2"])), RootLeafS)
RootleafS = intersect(as.character(unique(GeneFamilyDup[GeneFamilyDup$Var1 %in% RootLeafS, "Var2"])), RootS)
RootPollenLeafS = intersect(as.character(unique(GeneFamilyDup[GeneFamilyDup$Var1 %in% PollenLeafS, "Var2"])), RootS)
PollenLeafRootS = intersect(as.character(unique(GeneFamilyDup[GeneFamilyDup$Var1 %in% RootS, "Var2"])), PollenLeafS)
RootLeafPollenLeafS = intersect(as.character(unique(GeneFamilyDup[GeneFamilyDup$Var1 %in% PollenLeafS, "Var2"])), RootLeafS)
PollenLeafRootLeafS = intersect(as.character(unique(GeneFamilyDup[GeneFamilyDup$Var1 %in% RootLeafS, "Var2"])), PollenLeafS)


# Shared duplicated genes Large Scale Duplication (LSD) with no overlap
PollenDupLSDS = as.character(unique(paralogous[paralogous$Ref %in% PollenS,"Ref"]))
PollenLeafDupLSDS = as.character(unique(paralogous[paralogous$Ref %in% PollenLeafS,"Ref"]))
RootDupLSDS = as.character(unique(paralogous[paralogous$Ref %in% RootS,"Ref"]))
RootLeafDupLSDS = as.character(unique(paralogous[paralogous$Ref %in% RootLeafS,"Ref"]))
PollenPollenLSDS = intersect(as.character(unique(paralogous[paralogous$Ref %in% PollenS, "Anch"])), PollenS)
RootRootLSDS = intersect(as.character(unique(paralogous[paralogous$Ref %in% RootS, "Anch"])), RootS)
LeafLeafPollenLSDS = intersect(as.character(unique(paralogous[paralogous$Ref %in% PollenLeafS, "Anch"])), PollenLeafS)
LeafLeafRootLSDS = intersect(as.character(unique(paralogous[paralogous$Ref %in% RootLeafS, "Anch"])), RootLeafS)
RootPollenLSDS = intersect(as.character(unique(paralogous[paralogous$Ref %in% PollenS, "Anch"])), RootS)
PollenRootLSDS = intersect(as.character(unique(paralogous[paralogous$Ref %in% RootS, "Anch"])), PollenS)
LeafPollenLSDS = intersect(as.character(unique(paralogous[paralogous$Ref %in% PollenS, "Anch"])), PollenLeafS)
PollenRootLeafLSDS = intersect(as.character(unique(paralogous[paralogous$Ref %in% RootLeafS, "Anch"])), PollenS)
RootLeafPollenLSDS = intersect(as.character(unique(paralogous[paralogous$Ref %in% PollenS, "Anch"])), RootLeafS)
PollenleafLSDS = intersect(as.character(unique(paralogous[paralogous$Ref %in% PollenLeafS, "Anch"])), PollenS)
LeafRootLSDS = intersect(as.character(unique(paralogous[paralogous$Ref %in% RootS, "Anch"])), RootLeafS)
RootleafLSDS = intersect(as.character(unique(paralogous[paralogous$Ref %in% RootLeafS, "Anch"])), RootS)
RootPollenLeafLSDS = intersect(as.character(unique(paralogous[paralogous$Ref %in% PollenLeafS, "Anch"])), RootS)
PollenLeafRootLSDS = intersect(as.character(unique(paralogous[paralogous$Ref %in% RootS, "Anch"])), PollenLeafS)
RootLeafPollenLeafLSDS = intersect(as.character(unique(paralogous[paralogous$Ref %in% PollenLeafS, "Anch"])), RootLeafS)
PollenLeafRootLeafLSDS = intersect(as.character(unique(paralogous[paralogous$Ref %in% RootLeafS, "Anch"])), PollenLeafS)
#######################################################################################################################################################################################






#######################################################################################################################################################################################
# For entire module
# Also for Pollen-Root and Leaf-Leaf module senario since they also shared some APs
# Two heatmapes for each module pair of Pollen-Root and Leaf-Leaf, first one for all genes and the second one
# for genes which are not overpaled between two modules
#######################################################################################################################################################################################

Pollenleaf = RandomForest(classes = 2, data.exp.scale, PollenRoot, RootPollen, "Pollen", "Root", 50,
	title = "Pollen - Root modules overlapped APs", 
	address = "/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Plot/Anchorpoint/PollenRoot/PollenRootRF.pdf")

Pollenleaf = RandomForest(classes = 2, data.exp.scale, PollenRootS, RootPollenS, "Pollen", "Root", 50,
	title = "Pollen - Leaf modules Specific APs", 
	address = "/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Plot/Anchorpoint/PollenRoot/PollenRootRFS.pdf")



Leafleaf = RandomForest(classes = 2, data.exp.scale, PollenRootLeaf, RootPollenLeaf, "PollenLeaf", "RootLeaf", 50,
	title = "PollenLeaf - RootLeaf modules overlapped APs", 
	address = "/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Plot/Anchorpoint/PollenRoot/PollenRootLeafRF.pdf")

Leafleaf = RandomForest(classes = 2, data.exp.scale, PollenRootLeafS, RootPollenLeafS, "PollenLeaf", "RootLeaf", 50,
	title = "PollenLeaf - RootLeaf modules Specific APs", 
	address = "/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Plot/Anchorpoint/PollenRoot/PollenRootLeafRFS.pdf")
#######################################################################################################################################################################################




#######################################################################################################################################################################################
# Gene family and small-scale duplication as well as large-scale duplication for all four modules
#######################################################################################################################################################################################
GeneFamily = read.delim("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/GeneFamilyZhenFeb16TwoColumns.txt")
Pollen = as.character(read.table("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/Module/Pollen_90.txt")[,1])
PollenLeaf = as.character(read.table("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/Module/Leaf_90.txt")[,1])
Root = as.character(read.table("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/Module/AT1G50090_Root_90.txt")[,1])
RootLeaf = as.character(read.table("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/Module/AT3G19710_Leaf_90.txt")[,1])
ref = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/Reference_SD.txt", header = T)
rownames(ref) = ref[,1]
paralogous = read.table("~/Gene_duplication/Seed_genes/paralogous_NoFilter.txt", header = T)
paralogous$Ref = as.character(ref[paralogous$Ref,1])
paralogous$Anch = as.character(ref[paralogous$Anch,1])

GeneFamily$Duplication = "SSD"
for (i in 1:nrow(paralogous)) {
	GeneFamily[which((as.character(GeneFamily$Var1) == as.character(paralogous$Ref)[i] | as.character(GeneFamily$Var2) == as.character(paralogous$Ref)[i]) & (as.character(GeneFamily$Var1) == as.character(paralogous$Anch)[i] | as.character(GeneFamily$Var2) == as.character(paralogous$Anch)[i])), "Duplication"] = "LSD"
}


# Gene family size
GF = table(as.character(GeneFamily$GeneFamily))

# Gene family with atleast two genes
GFDup = names(GF)[GF > 1]
GeneFamilyDup = GeneFamily[GeneFamily$GeneFamily %in% GFDup,]
#######################################################################################################################################################################################





#######################################################################################################################################################################################
# Duplicated pairs (SSD-LSD) in four modules
#######################################################################################################################################################################################
DupRootLeaf = GeneFamilyDup[GeneFamilyDup$Var1 %in% RootLeaf & GeneFamilyDup$Var2 %in% Root, ]
DupRootLeaf = DupRootLeaf[!duplicated(DupRootLeaf[,1:2]),]
colnames(DupRootLeaf)[1:2] = c("Leaf", "Root")
DupRootLeaf$Region = "Overlap"
DupRootLeaf[DupRootLeaf$Root %in% Root & DupRootLeaf$Leaf %in% RootLeaf, "Region"] = "Specific"



length(unique(as.character(DupRootLeaf[DupRootLeaf$GeneFamily %in% "ORTHO000002_ORTHO001823_ORTHO013689", 1])))
length(unique(as.character(DupRootLeaf[DupRootLeaf$GeneFamily %in% "ORTHO000002_ORTHO001823_ORTHO013689", 2])))



DupPollenLeaf = GeneFamilyDup[GeneFamilyDup$Var1 %in% PollenLeaf & GeneFamilyDup$Var2 %in% Pollen, ]
DupPollenLeaf = DupPollenLeaf[!duplicated(DupPollenLeaf[,1:2]),]
colnames(DupPollenLeaf)[1:2] = c("Leaf", "Pollen")
DupPollenLeaf$Region = "Overlap"
DupPollenLeaf[DupPollenLeaf$Pollen %in% Pollen & DupPollenLeaf$Leaf %in% PollenLeaf, "Region"] = "Specific"


DupPollenRoot = GeneFamilyDup[GeneFamilyDup$Var1 %in% Root & GeneFamilyDup$Var2 %in% Pollen, ]
DupPollenRoot = DupPollenRoot[!duplicated(DupPollenRoot[,1:2]),]
colnames(DupPollenRoot)[1:2] = c("Root", "Pollen")
DupPollenRoot$Region = "Overlap"
DupPollenRoot[DupPollenRoot$Root %in% RootS & DupPollenRoot$Pollen %in% PollenS, "Region"] = "Specific"


DupLeafLeaf = GeneFamilyDup[GeneFamilyDup$Var1 %in% RootLeaf & GeneFamilyDup$Var2 %in% PollenLeaf, ]
DupLeafLeaf = DupLeafLeaf[!duplicated(DupLeafLeaf[,1:2]),]
colnames(DupLeafLeaf)[1:2] = c("RootLeaf", "PollenLeaf")
DupLeafLeaf$Region = "Overlap"
DupLeafLeaf[DupLeafLeaf$RootLeaf %in% RootLeafS & DupLeafLeaf$PollenLeaf %in% PollenLeafS, "Region"] = "Specific"



DupOverlap = GeneFamilyDup[GeneFamilyDup$Var1 %in% setdiff(Root, RootS) & GeneFamilyDup$Var2 %in% setdiff(PollenLeaf, PollenLeafS), ]
DupOverlap = DupOverlap[!duplicated(DupOverlap[,1:2]),]
colnames(DupOverlap)[1:2] = c("PollenRoot", "LeafLeaf")


DupMutwil = GeneFamilyDup[GeneFamilyDup$Var1 %in% mult.pollen & GeneFamilyDup$Var2 %in% mult.pollen, ]
DupMutwil = DupMutwil[!duplicated(DupMutwil[,1:2]),]
colnames(DupOverlap)[1:2] = c("PollenRoot", "LeafLeaf")



RandomForest = function(data, x, y, namex, namey, Top, address, classes, title, Modex, Modey, overlx, overly) {
	AnnotColor = c("Pollen" = "darkblue", "Root" = "red", "Leaf" = "chartreuse2", "Both" = "seashell4", "LSD" = "burlywood4", "SSD" = "gold", "Overlap" = "darkslategray1", "Specific" = "dodgerblue4", "PollenLeaf" = "darkgreen", "RootLeaf" = "darkorchid1", "PollenRoot" = "gray12", "LeafLeaf" = "green3")
	if(classes == 2) Data = data.frame(rbind.data.frame(data[x,], data[y,]), "class" = c(rep(namex, length(x)), rep(namey, length(y))))
	if(classes > 2) Data = data
	Result = randomForest(factor(class) ~ ., data = Data, ntrees = 20000)
	if(classes == 2) Var = data.frame(colnames(data), varUsed(Result))
	if(classes > 2) Var = data.frame(colnames(data[, -ncol(data)]), varUsed(Result))
	Var = Var[order(-Var[,2]),]
	
	if(is.null(overlx) & is.null(overly)){
		 annotationRow = data.frame("Module" = Data$class, "Mode" = c(Modex, Modey))
		 AnnotColorSelect = list("Module" = AnnotColor[names(AnnotColor) %in% names(table(Data$class))], "Mode" = AnnotColor[names(AnnotColor) %in% c(Modex, Modey)])
	}
	if(!(is.null(overlx) & is.null(overly))){
		annotationRow = data.frame("Module" = Data$class, "Mode" = c(Modex, Modey), "Overlap" = c(overlx, overly))
		AnnotColorSelect = list("Module" = AnnotColor[names(AnnotColor) %in% names(table(Data$class))], "Mode" = AnnotColor[names(AnnotColor) %in% c(Modex, Modey)],
		"Overlap" = AnnotColor[names(AnnotColor) %in% c(overlx, overly)])
	}

	ColorRamp = rgb(seq(0,1,length=256),seq(0,1,length=256),seq(1,0,length=256))
	ColorLevels = seq(-max(abs(Data[,-ncol(Data)])), max(abs(Data[,-ncol(Data)])), length=length(ColorRamp)+1)
	rownames(annotationRow) = rownames(Data)
	selected = Data[,as.numeric(rownames(Var[1:Top,]))]
	Col.n = GSM = c()
	for(i in 1:Top) {
		col.n = COR.annot.PO[names(COR.annot.PO) %in% Annot(colnames(selected)[i])][[1]]
		Col.n[i] = ifelse(length(col.n) == 1, col.n,
			ifelse(length(col.n) == 2, paste(col.n[1], col.n[2], sep = " | "),
				ifelse(length(col.n) == 3, paste(col.n[1], col.n[2], col.n[3], sep = " | "), colnames(selected)[i])))
		GSM[i] = COR.annot.GSM[names(COR.annot.GSM) %in% Annot(colnames(selected)[i])][[1]]
	}
	# RR = data.frame(Col.n, "A" = unlist(COR.annot.GSM[names(COR.annot.GSM) %in% Annot(colnames(selected))]))
	RR = data.frame(Col.n, GSM)
	RR = paste(RR[,1], RR[,2], sep = " * ")
	colnames(selected) = RR
	pheatmap(selected, cluster_row = F, clustering_method="ward.D2", annotation_row = annotationRow, border_color=NA, cluster_col = T,
		color = ColorRamp, breaks = ColorLevels, treeheight_col = 0, filename = address, fontsize = 8, main = title, annotation_colors = AnnotColorSelect,
		cellwidth = 2.5, cellheight = 6, fontsize_col = 2.5, fontsize_row = 6,  gaps_row = cumsum(table(factor(annotationRow[,1], levels = as.character(unique(annotationRow[,1]))))))
	# return(selected)
}


ModePollen = DupPollenLeaf[!duplicated(DupPollenLeaf$Pollen), "Duplication"]
ModeLeaf = DupPollenLeaf[!duplicated(DupPollenLeaf$Leaf), "Duplication"]
RandomForest(data = data.exp.scale, x = as.character(DupPollenLeaf[!duplicated(DupPollenLeaf$Pollen), "Pollen"]),
y = as.character(DupPollenLeaf[!duplicated(DupPollenLeaf$Leaf), "Leaf"]), "Pollen", "Leaf", title = "Pollen - Leaf module", classes = 2, Top = 30,
Modex = ModePollen, Modey = ModeLeaf, overlx = NULL, overly = NULL,
address = "/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Plot/Anchorpoint/RFTop10/TopConditionsPollenRFNew.pdf")


ModeRoot = DupRootLeaf[!duplicated(DupRootLeaf$Root), "Duplication"]
ModeLeaf = DupRootLeaf[!duplicated(DupRootLeaf$Leaf), "Duplication"]
RandomForest(data = data.exp.scale, x = as.character(DupRootLeaf[!duplicated(DupRootLeaf$Root), "Root"]),
y = as.character(DupRootLeaf[!duplicated(DupRootLeaf$Leaf), "Leaf"]), "Root", "Leaf", title = "Root - Leaf module", classes = 2, Top = 30,
Modex = ModeRoot, Modey = ModeLeaf, overlx = NULL, overly = NULL,
address = "/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Plot/Anchorpoint/RFTop10/TopConditionsRootRFNew.pdf")



ModePollenRoot = DupOverlap[!duplicated(DupOverlap$PollenRoot), "Duplication"]
ModeLeafLeaf = DupOverlap[!duplicated(DupOverlap$LeafLeaf), "Duplication"]
RandomForest(data = data.exp.scale, x = as.character(DupOverlap[!duplicated(DupOverlap$PollenRoot), "PollenRoot"]),
y = as.character(DupOverlap[!duplicated(DupOverlap$LeafLeaf), "LeafLeaf"]), "PollenRoot", "LeafLeaf", title = "PollenRoot - LeafLeaf module", classes = 2, Top = 30,
Modex = ModePollenRoot, Modey = ModeLeafLeaf, overlx = NULL, overly = NULL,
address = "/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Plot/Anchorpoint/RFTop10/TopConditionsOverlapRFNew.pdf")




ModeRoot = DupPollenRoot[!duplicated(DupPollenRoot$Root), "Duplication"]
ModePollen = DupPollenRoot[!duplicated(DupPollenRoot$Pollen), "Duplication"]
overlRoot = DupPollenRoot[!duplicated(DupPollenRoot$Root), "Region"]
overlPollen = DupPollenRoot[!duplicated(DupPollenRoot$Pollen), "Region"]
RandomForest(data = data.exp.scale, x = as.character(DupPollenRoot[!duplicated(DupPollenRoot$Root), "Root"]),
y = as.character(DupPollenRoot[!duplicated(DupPollenRoot$Pollen), "Pollen"]), "Root", "Pollen", title = "Root - Pollen module", classes = 2, Top = 30,
Modex = ModeRoot, Modey = ModePollen, overlx = overlRoot, overly = overlPollen,
address = "/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Plot/Anchorpoint/RFTop10/TopConditionsRootPollenRFNew.pdf")


ModeRootLef = DupPollenRoot[!duplicated(DupPollenRoot$Root), "Duplication"]
ModePollenLeaf = DupPollenRoot[!duplicated(DupPollenRoot$Pollen), "Duplication"]
overlRootLeaf = DupPollenRoot[!duplicated(DupPollenRoot$Root), "Region"]
overlPollenLeaf = DupPollenRoot[!duplicated(DupPollenRoot$Pollen), "Region"]

RandomForest(data = data.exp.scale, x = as.character(DupPollenRoot[!duplicated(DupPollenRoot$Root), "Root"]),
y = as.character(DupPollenRoot[!duplicated(DupPollenRoot$Pollen), "Pollen"]), "Root", "Pollen", title = "Root - Pollen module", classes = 2, Top = 30,
Modex = ModeRoot, Modey = ModePollen, overlx = overlRoot, overly = overlPollen,
address = "/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Plot/Anchorpoint/RFTop10/TopConditionsRootPollenRFNew.pdf")




ModeRootLeaf = DupLeafLeaf[!duplicated(DupLeafLeaf$RootLeaf), "Duplication"]
ModePollenLeaf = DupLeafLeaf[!duplicated(DupLeafLeaf$PollenLeaf), "Duplication"]
overlRootLeaf = DupLeafLeaf[!duplicated(DupLeafLeaf$RootLeaf), "Region"]
overlPollenLeaf = DupLeafLeaf[!duplicated(DupLeafLeaf$PollenLeaf), "Region"]

RandomForest(data = data.exp.scale, x = as.character(DupLeafLeaf[!duplicated(DupLeafLeaf$RootLeaf), "RootLeaf"]),
y = as.character(DupLeafLeaf[!duplicated(DupLeafLeaf$PollenLeaf), "PollenLeaf"]), "RootLeaf", "PollenLeaf", title = "RootLeaf - PollenLeaf module", classes = 2, Top = 30,
Modex = ModeRootLeaf, Modey = ModePollenLeaf, overlx = overlRootLeaf, overly = overlPollenLeaf,
address = "/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Plot/Anchorpoint/RFTop10/TopConditionsLeafLeafRFNew.pdf")
#######################################################################################################################################################################################





Randomization = function(x, y) {
	Sample = as.character(ref$Gene[sample(1:nrow(ref), c(length(x) + length(y)))])
	SampleRoot = Sample[1:length(x)]
	SampleRootLeaf = Sample[c(length(x) + 1):c(length(x) + length(y))]
	DupSampleRootLeaf = GeneFamilyDup[GeneFamilyDup$Var1 %in% SampleRootLeaf & GeneFamilyDup$Var2 %in% SampleRoot, ]
	DupSampleRootLeaf = DupSampleRootLeaf[!duplicated(DupSampleRootLeaf[,1:2]),]
	Result = table(DupSampleRootLeaf$Duplication)
	return(Result)
}


RR.RootLeaf = foreach(i =1:1000, .combine = rbind) %do% {
	print(i)
	Randomization(Root, RootLeaf)
}
