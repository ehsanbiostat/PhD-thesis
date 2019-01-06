

#################################################################
#
# Deduce function of APs based on highly similar genes, gene score correlation, and GSEA
#
#################################################################

library(igraph)
library(foreach)
library(doMC)
options(width = 200)
# registerDoMC(2)
k = 6
condition.names = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/condition_name.txt", header = T)
NAME = names(table(condition.names))
NAME.spe = NAME[k]
source("/ngsprojects/iwt/ehsab/Gene_duplication/Scripts/General/myImagePlot.R")
source("/ngsprojects/iwt/ehsab/Gene_duplication/Scripts/General/myImagePlot_sym.R")
ref = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/Reference_SD.txt", header = T)
rownames(ref) = ref[,1]


LeafRoot = read.table("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/Module/AT3G19710_Leaf_90.txt")
Root = read.table("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/Module/AT1G50090_Root_90.txt")
AnchRoot = read.delim("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/Anchorpoint/AT1G50090_AT3G19710_52_90.txt")



Leaf = read.table("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/Module/Leaf_90.txt")
Pollen = read.table("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/Module/Pollen_90.txt")
Anch = read.delim("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/Anchorpoint/Geneset1_2_53_90.txt.txt")

PollenRef = ref[as.character(Pollen[,1]), 2]
PollenAnch = unique(Anch$Pollen)

LeafRef = ref[as.character(Leaf[,1]), 2]
LeafAnch = unique(Anch$Leaf)


RootRef = ref[as.character(Root[,1]), 2]
RootAnch = unique(AnchRoot$Root)


LeafRootRef = ref[as.character(LeafRoot[,1]), 2]
LeafRootAnch = unique(AnchRoot$Leaf)


# CondPollen = matrix(NA, length(PollenRef), sum(condition.names == NAME.spe))
# CondLeaf = matrix(NA, length(LeafRef), sum(condition.names == NAME.spe))

# j = 0
# for(i in PollenRef) {
# 	j = j + 1
# 	load(paste("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Condition_score/Compendium_wise/development/Raw/C",PollenRef[j],".RData", sep = ""))
# 	CondPollen[j,] = cond_scores[,1]
# }

# rownames(CondPollen) = as.character(Pollen[,1])


# PollenCorr = foreach(i = 1:length(PollenAnch), .combine = cbind) %do% {
# 	AA = foreach(j = 1:nrow(Pollen), .combine = c) %do% {
# 		cor(CondPollen[as.character(PollenAnch[i]),], CondPollen[j,])
# 	}
# 	AA
# }

# colnames(PollenCorr) = as.character(PollenAnch)
# rownames(PollenCorr) = as.character(Pollen[,1])

# PollenCorr = as.data.frame(PollenCorr)
# PollenCorr["Gene"] = rownames(PollenCorr)




# j = 0
# for(i in LeafRef) {
# 	j = j + 1
# 	load(paste("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Condition_score/Compendium_wise/development/Raw/C",LeafRef[j],".RData", sep = ""))
# 	CondLeaf[j,] = cond_scores[,1]
# }

# rownames(CondLeaf) = as.character(Leaf[,1])


# LeafCorr = foreach(i = 1:length(LeafAnch), .combine = cbind) %do% {
# 	AA = foreach(j = 1:nrow(Leaf), .combine = c) %do% {
# 		cor(CondLeaf[as.character(LeafAnch[i]),], CondLeaf[j,])
# 	}
# 	AA
# }

# colnames(LeafCorr) = as.character(LeafAnch)
# rownames(LeafCorr) = as.character(Leaf[,1])

# LeafCorr = as.data.frame(LeafCorr)
# LeafCorr["Gene"] = rownames(LeafCorr)




##################################
# Gene score as GSEA scores
##################################

AnchPollenRef = unique(Anch$Ref)
AnchLeafRef = unique(Anch$References)
AnchRootRef = unique(AnchRoot$Ref)
AnchLeafRootRef = unique(AnchRoot$References)


GeneScore = function(Ref, AnchRef, Anch, module) { # Generate gene score for anchorpoints and genes within the module
	GeneScore = matrix(NA, 19285, length(AnchRef))
	j = 0
	for(i in AnchRef) {
		j = j + 1
		load(paste("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Gene_score/Compendium_wise/development/Raw/G",AnchRef[j],".RData", sep = ""))
		GeneScore[,j] = gene_scores[,1]
	}
	colnames(GeneScore) = as.character(Anch)
	rownames(GeneScore) = as.character(ref[,1])
	GeneScore = as.data.frame(GeneScore)
	GeneScore["Gene"] = rownames(GeneScore)
	return(GeneScore)
}

GeneScoreRoot = GeneScore(RootRef, AnchRootRef, RootAnch, Root)
GeneScoreLeafRoot = GeneScore(LeafRootRef, AnchLeafRootRef, LeafRootAnch, LeafRoot)

GeneScorePollen = GeneScore(PollenRef, AnchPollenRef, PollenAnch, Pollen)
GeneScoreLeaf = GeneScore(LeafRef, AnchLeafRef, LeafAnch, Leaf)

# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# GO analysis
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
library(igraph)
library(foreach)
library(doMC)
# registerDoMC(2)
library(ape)
options(width = 200)

go = read.delim("/ngsprojects/iwt/ehsab/Prioritization/prioritization/Combined_Network/New/GO_OCT2015.txt", header = F)
clas = c("C", "F", "P")
k = c(2,3)
go = go[go[, 4] %in% clas[k], ]
go.limited = go[!(go$V2 %in%  rownames(data.frame(sort(table(as.matrix(go$V2)), decreasing=T)[1:20]))),]
goMap = go[!(duplicated(go[,c(2:3)])),]
# go = go[go[, 4] == clas[k] & (go[, 5] %in% c("EXP", "IDA","IPI", "IMP", "IGI", "IEP")), ] # only experimental evidence codes





# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# GSEA analysis
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
library(piano)
# myGsc = loadGSC(as.matrix(go.limited[, c("V1", "V2")]))
registerDoMC(20)


GSEA = function(Scores, method, FDR, module, TOP) { #Performs GSEA 
	T = foreach (i = 1:c(ncol(Scores) - 1)) %dopar% {
		print(i)
		if(is.na(TOP)) {
			myStats = Scores[,i]
			names(myStats) = as.character(rownames(Scores))
		} else {
			select = which(rank(-Scores[,i]) <= TOP)
			GOSelect = unique(as.character(go.limited[go.limited$V1 %in% rownames(Scores)[select],"V3"]))
			go.limitedSelect = go.limited[as.character(go.limited$V3) %in% GOSelect,]
			GeneSelect = unique(as.character(go.limitedSelect[,1]))
			myStats = Scores[rownames(Scores) %in% GeneSelect,i]
			names(myStats) = as.character(rownames(Scores)[rownames(Scores) %in% GeneSelect])
			myGsc = loadGSC(as.matrix(go.limitedSelect[, c("V1", "V2")]))
			# select = which(Scores[,i] > 0)
			# myStats = 1/(1-Scores[select, i])
			# names(myStats) = as.character(rownames(Scores)[select])
		}
		gsaRes <- runGSA(myStats, gsc=myGsc, geneSetStat = method, gsSizeLim=c(5,500), ncpus=1)
		gsaRes
	}
	
	GSAtable = lapply(T, GSAsummaryTable)
	
	final = foreach (i = 1:length(GSAtable), .combine = rbind) %dopar% {
		print(i)
		x = GSAtable[[i]]
		x = x[x[, "p adj (non-dir.)"] < FDR,]
		if(nrow(x) > 0) x["Gene"] = colnames(Scores)[i]
		if(nrow(x) > 0) x[, c("Name", "Genes (tot)", "p adj (non-dir.)", "Gene")]
	}
	final["module"] = module
	return(final)
}

GSEAPollen = GSEA(GeneScorePollen	, "mean", 0.1, "Pollen", TOP = 50)
GSEAPollenLeaf = GSEA(GeneScoreLeaf, "mean", 0.1, "Leaf", TOP = 50)

GSEARoot = GSEA(GeneScoreRoot, "mean", 0.1, "Root", TOP = 50)
GSEARootLeaf = GSEA(GeneScoreLeafRoot, "mean", 0.1, "Leaf", TOP = 50)


write.table(rbind(GSEAPollen, GSEAPollenLeaf), file = "/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/Anchorpoint/GSEA_PollenLeaf_Top50_500.txt", row.names = F, quote = F , sep = "\t")
write.table(rbind(GSEARoot, GSEARootLeaf), file = "/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/Anchorpoint/GSEA_RootLeaf_Top50_500.txt", row.names = F, quote = F , sep = "\t")




#################################################################################################################################################################################
#  Heatmap and clustering on both row and column of GSEA matrix, columns are GO terms and rows are genes
#################################################################################################################################################################################

library(gplots)
library(made4)
library(igraph)
library(NMF)
nmf.options(grid.patch=TRUE)
# install.packages("pheatmap", repos = "http://lib.ugent.be/CRAN/")
library(pheatmap)
#1E3C5C
PollenLeaf = read.delim("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/Anchorpoint/GSEA_PollenLeaf_Top50_500.txt")
RootLeaf = read.delim("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/Anchorpoint/GSEA_RootLeaf_Top50_500.txt")
RootLeafGO = merge(goMap, RootLeaf, by.x = "V2", by.y = "Name")[,c(1,3,6,7,8,9)]
RootLeafGO = RootLeafGO[order(RootLeafGO$module),]
PollenLeafGO = merge(goMap, PollenLeaf, by.x = "V2", by.y = "Name")[,c(1,3,6,7,8,9)]
PollenLeafGO = PollenLeafGO[order(PollenLeafGO$module),]
RootLeafGO["module"] = ifelse(RootLeafGO$module == "Leaf", "Root-Leaf", "Root")
PollenLeafGO["module"] = ifelse(PollenLeafGO$module == "Leaf", "Pollen-Leaf", "Pollen")
PollenRootLeafGO = rbind(PollenLeafGO, RootLeafGO)
PollenRootLeaf = rbind(PollenLeaf, RootLeaf)

PollenLeaf[PollenLeaf$Gene %in% unique(intersect(as.character(PollenLeaf$Gene), as.character(RootLeaf$Gene)))[i],]
RootLeaf[RootLeaf$Gene %in% unique(intersect(as.character(PollenLeaf$Gene), as.character(RootLeaf$Gene)))[i],]

# Function to convert GSEA output to a corresponding matrix in order to visualize it by heatmap and hirerchical clustering in the next setps
Mydata = function(x, FDR) {
	x = x[x$p.adj..non.dir.. < FDR,]
	k = 0
	if (ncol(x) == 6) k = 1
	g = graph.edgelist(as.matrix(x[, c(1+k, 4+k)]), directed=TRUE)
	E(g)$weight = -log(x[, "p.adj..non.dir.."] + 0.00001)
	AJ = get.adjacency(g, attr="weight", sparse=F)
	Aj1 = AJ[,-c(which(colnames(AJ) %in% x[,1+k]))]
	mydata = Aj1[-c(which(rownames(Aj1) %in% x[,4+k])),]
	return(mydata)
}



ReorderMatrix = function(x) { 
	# set the custom distance and clustering functions, per your example
	hclustfunc <- function(x) hclust(x, method="ward.D2")
	distfunc <- function(x) dist(x, method="euclidean")

	# perform clustering on rows and columns
	cl.row <- hclustfunc(distfunc(t(mydata)))
	cl.col <- hclustfunc(distfunc(mydata))
	mydata = mydata[cl.col$order, cl.row$order]
	return(mydata)
}

# extract cluster assignments; i.e. k=8 (rows) k=5 (columns)
gr.row <- cutree(cl.row, 8)
gr.col <- cutree(cl.col, 8)

# require(RColorBrewer)
col1 <- brewer.pal(8, "Dark2")
col2 <- brewer.pal(8, "Dark2")


# jpeg("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/GSEA_RootLeaf.jpeg", width = 13, height = 7, units = "in", res = 800, pointsize=8)
# heatmap.2(t(mydata), hclustfun=hclustfunc, distfun=distfunc,
# 	RowSideColors=col1[gr.row], ColSideColors=col2[gr.col], trace = "none", scale = NULL, keysize=0.7, srtCol=45, adjCol = c(1,0), col= colorpanel(5, "white", "black"),
# 	sepcolor='yellow', sepwidth=0.05)
# dev.off()


Type = "Top50_500"

mydata = Mydata(RootLeaf, 0.1)
write.table(ReorderMatrix(mydata), file =paste("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/GO/GSEA/GSEA_RootLeaf",Type,".txt", sep = ""), quote = F, sep = "\t")
annotation = data.frame("Module" = ifelse(colnames(mydata) %in% AnchRoot$Root, "Root", "Leaf"), "Noverlap" = ifelse(colnames(mydata) %in% AnchRoot$Root & !(colnames(mydata) %in% Anch$Pollen), "Root",
	ifelse(colnames(mydata) %in% AnchRoot$Root & colnames(mydata) %in% Anch$Pollen, "Root-Pollen",
		ifelse(colnames(mydata) %in% AnchRoot$Leaf & colnames(mydata) %in% Anch$Leaf, "LeafLeaf", "Leaf"))))
rownames(annotation) = colnames(mydata)
# pdf("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/GSEA_RootLeaf.pdf", width = 20, height = 15, pointsize=12)
# aheatmap(t(mydata), color = c("white","black"), Rowv = TRUE, Colv = TRUE, labCol = NULL, distfun = "euclidean", hclustfun = "ward",
# 	annRow = annotation, annColors = list(c("blue", "red")), cellwidth = 3, cellheight = 4, legend = TRUE)
# dev.off()



pheatmap(t(mydata), color = colorRampPalette(c("white", "black"))(10), clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean", clustering_method = "ward.D2",
	labCol = NULL, annotation_row = annotation, cellwidth = 3, cellheight = 4, legend = TRUE, fontsize_col = 3, fontsize_row = 4,
	border_color = NA, cutree_rows = 4, cutree_cols = 4, , annotation_colors = list(Module = c(Root = "darkblue", Leaf = "red"), Noverlap = c(Root = "darkblue", Leaf = "red", LeafLeaf = "green", "Root-Pollen" = "yellow")),
	filename = paste("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/GO/GSEA/GSEA_RootLeaf",Type,".pdf", sep = ""))


mydata = Mydata(PollenLeaf, 0.1)
write.table(ReorderMatrix(mydata), file =paste("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/GO/GSEA/GSEA_PollenLeaf",Type,".txt", sep = ""), quote = F, sep = "\t")
annotation = data.frame("Module" = ifelse(colnames(mydata) %in% Anch$Pollen, "Pollen", "Leaf"), "Noverlap" = ifelse(!(colnames(mydata) %in% AnchRoot$Root) & colnames(mydata) %in% Anch$Pollen, "Pollen",
	ifelse(colnames(mydata) %in% AnchRoot$Root & colnames(mydata) %in% Anch$Pollen, "Root-Pollen",
		ifelse(colnames(mydata) %in% AnchRoot$Leaf & colnames(mydata) %in% Anch$Leaf, "LeafLeaf", "Leaf"))))
rownames(annotation) = colnames(mydata)
pheatmap(t(mydata), color = colorRampPalette(c("white", "black"))(10), clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean", clustering_method = "ward.D2",
	labCol = NULL, annotation_row = annotation, annotation_colors = list(Module = c(Pollen = "darkblue", Leaf = "red"), Noverlap = c(Pollen = "darkblue", Leaf = "red", LeafLeaf = "green", "Root-Pollen" = "yellow")), cellwidth = 3, cellheight = 4, legend = TRUE, fontsize_col = 3, fontsize_row = 4,
	border_color = NA, cutree_rows = 5, cutree_cols = 5,
	filename = paste("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/GO/GSEA/GSEA_PollenLeaf",Type,".pdf", sep = ""))


mydata = Mydata(PollenRootLeafGO, 0.1)
write.table(ReorderMatrix(mydata), file =paste("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/GO/GSEA/GSEA_PollenRootLeafGO",Type,".txt", sep = ""), quote = F, sep = "\t")
annotation = data.frame("Module" = factor(ifelse(colnames(mydata) %in% Anch$Pollen, "Pollen", ifelse(colnames(mydata) %in% AnchRoot$Root, "Root", ifelse(colnames(mydata) %in% AnchRoot$Leaf, "Root-Leaf", "Pollen-Leaf")))))
rownames(annotation) = colnames(mydata)
pheatmap(t(mydata), color = colorRampPalette(c("white", "black"))(10), clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean", clustering_method = "ward.D2",
	labCol = NULL, annotation_row = annotation, annotation_colors = list(Module = c(Pollen = "darkblue", Root = "red", "Root-Leaf" = "green", "Pollen-Leaf" ="yellow")), cellwidth = 3, cellheight = 4, legend = TRUE, fontsize_col = 3, fontsize_row = 4,
	border_color = NA, cutree_rows = 6, cutree_cols = 6,
	filename = paste("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/GO/GSEA/GSEA_PollenRootLeafGO",Type,".pdf", sep = ""))


mydata = Mydata(PollenRootLeaf, 0.05)
write.table(ReorderMatrix(mydata), file =paste("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/GO/GSEA/GSEA_PollenRootLeaf",Type,".txt", sep = ""), quote = F, sep = "\t")
annotation = data.frame("Module" = ifelse(colnames(mydata) %in% Anch$Pollen, "Pollen", ifelse(colnames(mydata) %in% AnchRoot$Root, "Root", ifelse(colnames(mydata) %in% AnchRoot$Leaf, "Root-Leaf", "Pollen-Leaf"))))
rownames(annotation) = colnames(mydata)
pheatmap(t(mydata), color = colorRampPalette(c("white", "black"))(10), clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean", clustering_method = "ward.D2",
	labCol = NULL, annotation_row = annotation, annotation_colors = list(Module = c(Pollen = "darkblue", Root = "red", "Root-Leaf" = "green", "Pollen-Leaf" ="yellow")), cellwidth = 3, cellheight = 4, legend = TRUE, fontsize_col = 3, fontsize_row = 4,
	border_color = NA, cutree_rows = 6, cutree_cols = 7, annotation_col = annotationCol,
	filename = paste("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/GO/GSEA/GSEA_PollenRootLeaf",Type,".pdf", sep = ""))



G1 = rownames(mydata)[1:14]
G2 = rownames(mydata)[15:16]
G3 = rownames(mydata)[17:24]
G4 = rownames(mydata)[25:55]
G5 = rownames(mydata)[602:619]
G6 = rownames(mydata)[590:601]

GG = rbind(data.frame("GOTerms" = rownames(mydata)[1:14], "Group" = "G1"),
	data.frame("GOTerms" = rownames(mydata)[15:16], "Group" = "G2"),
	data.frame("GOTerms" = rownames(mydata)[17:24], "Group" = "G3"),
	data.frame("GOTerms" = rownames(mydata)[25:55], "Group" = "G4"),
	data.frame("GOTerms" = rownames(mydata)[602:619], "Group" = "G5"),
	data.frame("GOTerms" = rownames(mydata)[590:601], "Group" = "G6")
	)

UU = rbind(GG, data.frame("GOTerms" = setdiff(rownames(mydata), as.character(GG$GOTerms)), "Group" = "Other"))
annotationCol = data.frame(UU[,2])
rownames(annotationCol) = as.character(UU$GOTerms)
colnames(annotationCol) = "GO Groups"


GOO = matrix(0, ncol(mydata), length(table(GG$Group)))
rownames(GOO) = colnames(mydata)
colnames(GOO) = names(table(GG$Group))


for(i in 1:ncol(mydata)) {
	GOO[colnames(mydata)[i],unique(as.character(GG$Group[GG$GOTerms %in% rownames(mydata)[mydata[,i] > 0]]))] = 1
}

GG1 = data.frame("G1" = 0, "G2" = 0, "G3" = 0, "G4" = 0, "G5" = 0, "G6" = 0)
rownames(GG1) = "AT3G59440"
GOO = rbind(GOO, GG1)
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
	rownames(annotationRow) = rownames(Data)
	annotationRow = cbind(annotationRow, rbind(GOO[x,], GOO[y,]))
	ColorRamp = rgb(seq(0,1,length=256),seq(0,1,length=256),seq(1,0,length=256))
	ColorLevels = seq(-max(abs(Data[,-ncol(Data)])), max(abs(Data[,-ncol(Data)])), length=length(ColorRamp)+1)
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
		color = ColorRamp, breaks = ColorLevels, treeheight_col = 0,  fontsize = 5, main = title, filename = address,
		cellwidth = 6, cellheight = 5, fontsize_col = 5, fontsize_row = 5,  gaps_row = cumsum(table(factor(annotationRow[,1], levels = as.character(unique(annotationRow[,1]))))))
	return(selected)
}




PollenLeaf = RandomForest(classes = 2, data.exp.scale, as.character(AnchPollen$Pollen), as.character(AnchPollen$Leaf), "Pollen", "Leaf", 50,
	title = "Pollen - Leaf modules", 
	address = "/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Plot/Anchorpoint/RFTop10/TopConditionsPollenRF1.pdf")

RootLeaf = RandomForest(classes = 2, data.exp.scale, as.character(AnchRoot$Root), as.character(AnchRoot$Leaf), "Root", "Leaf", 50,
	title = "Root - Leaf modules", 
	address = "/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Plot/Anchorpoint/RFTop10/TopConditionsRootRF1.pdf")









# Extracting clusters both for anchotpoints and GO enriched terms
# For Root-leaf module
# Anch.clusters = as.data.frame(gr.row)
# Anch.clusters["gr.row"] = ifelse(gr.row == 7 , "brown", ifelse(gr.row == 6 , "yellow", ifelse(gr.row == 5, "lightgreen", ifelse(gr.row == 4, "Pink", ifelse(gr.row == 3, "purple", ifelse(gr.row == 2, "orange", ifelse(gr.row == 1, "darkgreen", "gray")))))))
# Anch.clusters["gr.row"] = ifelse(gr.row == 5, "lightgreen", ifelse(gr.row == 4, "Pink", ifelse(gr.row == 3, "purple", ifelse(gr.row == 2, "orange", "darkgreen"))))
# GSEA = data.frame("Cluster" = ifelse(gr.col == 1, "darkgreen", ifelse(gr.col == 2, "orange", ifelse(gr.col == 3, "purple", ifelse(gr.col == 4, "pink", ifelse(gr.col == 5, "lightgreen", ifelse(gr.col == 6 , "yellow", ifelse(gr.col == 7 , "brown", "gray"))))))))
# Anch.clusters["module"] = "Root"
# Anch.clusters["module"] = "Pollen"
# Anch.clusters[rownames(Anch.clusters) %in% as.character(Anch$Leaf), "module"] = "Leaf"
# Anch.clusters["Anchorpoint"] = rownames(Anch.clusters)
# table(Anch.clusters$module, Anch.clusters$gr.row)
# write.table(Anch.clusters, file = "/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/GO/GSEA/APsRootLeaf.txt", row.names = F, quote = F, sep = "\t")
# write.table(GSEA, file = "/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/GO/GSEA/GSEARootLeaf.txt", row.names = T, quote = F, sep = "\t")






goExp = go[go$V5 %in% c("IDA","IEP","IPI","IMP","IGI"),]


PollenGene = as.character(unique(PollenLeaf$Gene))


pollenleaf = foreach (i = 1:length(PollenGene), .combine = rbind) %do% {
	x = as.character(goExp[goExp$V1 %in% PollenGene[i], "V2"])
	y = PollenLeaf[PollenLeaf$Gene %in% PollenGene[i], "Name"]
	Int = intersect(x, y)
	GOterm = "NA"
	if(length(x) > 0 & length(Int) == 0) GOterm = paste(x, collapse=" | ")
	if(length(x) > 0 & length(Int) > 0) GOterm = paste(c(paste(Int, "*", sep = ""), setdiff(x, Int)), collapse= " | ")
	data.frame("Gene" = PollenGene[i], "EXP_GO" = length(x), "GSEA" = length(y),
		"Jaccard" = length(intersect(x, y))/length(union(x, y)), "Module" = as.character(PollenLeaf[PollenLeaf$Gene %in% PollenGene[i], "module"][1]), "EXPGOterm" = GOterm)
}

sum(pollenleaf[pollenleaf$Exp > 0, "Jaccard"] > 0)/length(pollenleaf[pollenleaf$Exp > 0, "Jaccard"])

write.table(pollenleaf, file = "/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/GO/GSEA/Overlap_GSEA_EXP_PollenLeaf.txt", row.names = F, sep = "\t", quote = F)


PollenGene = as.character(unique(RootLeaf$Gene))
rootleaf = foreach (i = 1:length(PollenGene), .combine = rbind) %do% {
	x = goExp[goExp$V1 %in% PollenGene[i], "V2"]
	y = RootLeaf[RootLeaf$Gene %in% PollenGene[i], "Name"]
	Int = intersect(x, y)
	GOterm = "NA"
	if(length(x) > 0 & length(Int) == 0) GOterm = paste(x, collapse=" | ")
	if(length(x) > 0 & length(Int) > 0) GOterm = paste(c(paste(Int, "*", sep = ""), setdiff(x, Int)), collapse= " | ")
	data.frame("Gene" = PollenGene[i], "Exp" = length(x), "GSEA" = length(y),
		"Jaccard" = length(intersect(x, y))/length(union(x, y)), "Module" = as.character(RootLeaf[RootLeaf$Gene %in% PollenGene[i], "module"][1]), "EXPGOterm" = GOterm)
}

sum(rootleaf[rootleaf$Exp > 0, "Jaccard"] > 0)/length(rootleaf[rootleaf$Exp > 0, "Jaccard"])
nrow(rootleaf[rootleaf$Exp > 0 & rootleaf$Module == "Root" & rootleaf$Jaccard > 0, ])


write.table(rootleaf, file = "/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/GO/GSEA/Overlap_GSEA_EXP_RootLeaf.txt", row.names = F, sep = "\t", quote = F)

