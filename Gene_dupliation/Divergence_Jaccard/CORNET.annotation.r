library(igraph)
library(foreach)
library(doMC)
registerDoMC(16)

## Read condition names
condition.names = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/condition_name.txt", header = T)
NAME = names(table(condition.names))
PO = read.delim("/group/biocomp/projects/CORNET/MA_annotation/cornet_allMA_desc_120509.txt")
k = 1
NAME.spe = NAME[k]

## Condition scores
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

## Load generated graph from SA result
load(paste("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Gene_score/Compendium_wise/",NAME[k],"/Graph_V2.RData", sep = ""))


Result = read.table(paste("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Local/Compendium_wise/",NAME[k],"/Score99/Scores_V2/Anch_Score.txt", sep = ""), header = T)
ref = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/Reference_SD.txt", header = T)
# paralogous = read.table("~/Gene_duplication/Seed_genes/paralogous_total.txt", header = T)
paralogous = read.table("~/Gene_duplication/Seed_genes/paralogous_NoFilter.txt", header = T)
PARA = unique(paralogous$Ref)
PARA.recent = paralogous[paralogous$WGDevent == "recent", 1]
PARA.recent = unique(PARA.recent)
PARA.ancient = paralogous[paralogous$WGDevent == "ancient", 1]
PARA.ancient = unique(PARA.ancient)
Result["Jaccard2"] = (Result$Cluster.overlap + 0.0001)/(Result$Gene.cluster.size + Result$Paralogous.cluster.size)


## Gene expression
data.exp = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/Integrated_SD.txt", colClasses = "numeric", nrows = 19285, comment.char = "", header = T)
data.exp = data.exp[, condition.names == NAME.spe]

Module.annot = matrix(NA, nrow(Result) * 20, 4)
colnames(Module.annot) = c("Annot", "Jaccard", "Direction", "Gene")
Module.annot = as.data.frame(Module.annot)
# Module.annot = foreach (j = 1:nrow(Result), .combine = rbind) %dopar% {

m = -19
for (j in 1:nrow(Result)) {	
	m = m + 20
 	print(j)
	memb.comun = neighborhood(g, 1, Result$Gene[j], mode = c("out"))[[1]]
	# memb.comun.para = neighborhood(g, 1, Result$Paralogous[j], mode = c("out"))[[1]]
	
	CondRank1 = foreach(i = 1:length(memb.comun), .combine = rbind) %do%{
		rank(total_cond_scores[, memb.comun[i]])
	}
	if (length(memb.comun) > 1) rownames(CondRank1) = c(1:length(memb.comun))
	CondRank = rbind(c(1:dim(total_cond_scores)[1]), CondRank1)
	colnames(CondRank) = colnames(data.exp)
	ORDER = 1
	CondRank0 = CondRank[,order(CondRank[ORDER+1, ], decreasing = T)]
	RR = rbind.data.frame(cbind.data.frame("Annot" = colnames(CondRank0)[1:10], "Jaccard" = Result[j, "Jaccard2"], "Direction" = "Up"),
		cbind.data.frame("Annot" = colnames(CondRank0)[c(nrow(total_cond_scores)-10):nrow(total_cond_scores)], "Jaccard" = Result[j, "Jaccard2"], "Direction" = "Down"))
	RR = cbind(RR, "Gene" = Result$Gene[j])
	Module.annot[m:c(m+20), 1] = as.character(RR[,1])
	Module.annot[m:c(m+20), 2] = as.data.frame(RR[,2])	
	Module.annot[m:c(m+20), 3] = as.character(RR[,3])
	Module.annot[m:c(m+20), 4] = as.data.frame(RR[,4])
}

# Remove corresponding anchorpoints which have less than three coexpressed genes in their module
Module.annot = Module.annot[!(Module.annot$Gene %in% Result[Result[,"Gene.cluster.size"] < 3, "Gene"]), ]


# Divide the results into two groups for up and down regulated conditions
Module.annot.up = Module.annot[Module.annot$Direction == "Up", ]
Module.annot.down = Module.annot[Module.annot$Direction == "Down", ]

Module.annot.up.recent = Module.annot.up[Module.annot.up$Gene %in% PARA.recent,]
Module.annot.down.recent = Module.annot.down[Module.annot.down$Gene %in% PARA.recent,]
Module.annot.up.ancient = Module.annot.up[Module.annot.up$Gene %in% PARA.ancient,]
Module.annot.down.ancient = Module.annot.down[Module.annot.down$Gene %in% PARA.ancient,]
# --------------------------------------------------------------------------------------------------------------------------------------------------
# Scores
# --------------------------------------------------------------------------------------------------------------------------------------------------
# Recent anchorpoints with up regulated conditions
Annot.score.up.recent.mean = aggregate(Module.annot.up.recent$Jaccard ~ Module.annot.up.recent$Annot, FUN = mean)
Annot.score.up.recent.med = aggregate(Module.annot.up.recent$Jaccard ~ Module.annot.up.recent$Annot, FUN = median)
Annot.score.up.recent.sd = aggregate(Module.annot.up.recent$Jaccard ~ Module.annot.up.recent$Annot, FUN = sd)
AA = table(Module.annot.up.recent$Annot)
AA = AA[names(AA) %in% Annot.score.up.recent.med[,1]]
Annot.score.up.recent = cbind.data.frame(Annot.score.up.recent.mean, "Freq" = AA, "median" = Annot.score.up.recent.med[,2], "SD" = Annot.score.up.recent.med[,2])
Annot.score.up.recent = Annot.score.up.recent[order(Annot.score.up.recent$median), ]
Annot.score.up.recent[, 1] = sub("X", "", Annot.score.up.recent[,1])
Annot.score.up.recent.3 = Annot.score.up.recent[Annot.score.up.recent$Freq > 3,]
Annot.score.up.recent.3["Rank"] = 1:nrow(Annot.score.up.recent.3) 

# Recent anchorpoints with down regulated conditions
Annot.score.down.recent.mean = aggregate(Module.annot.down.recent$Jaccard ~ Module.annot.down.recent$Annot, FUN = mean)
Annot.score.down.recent.med = aggregate(Module.annot.down.recent$Jaccard ~ Module.annot.down.recent$Annot, FUN = median)
Annot.score.down.recent.sd = aggregate(Module.annot.down.recent$Jaccard ~ Module.annot.down.recent$Annot, FUN = sd)
AA = table(Module.annot.down.recent$Annot)
AA = AA[names(AA) %in% Annot.score.down.recent.med[,1]]
Annot.score.down.recent = cbind.data.frame(Annot.score.down.recent.mean, "Freq" = AA, "median" = Annot.score.down.recent.med[,2], "SD" = Annot.score.down.recent.med[,2])
Annot.score.down.recent = Annot.score.down.recent[order(Annot.score.down.recent$median), ]
Annot.score.down.recent[, 1] = sub("X", "", Annot.score.down.recent[,1])
Annot.score.down.recent.3 = Annot.score.down.recent[Annot.score.down.recent$Freq > 3,]
Annot.score.down.recent.3["Rank"] = 1:nrow(Annot.score.down.recent.3) 

# Ancient anchorpoints with up regulated conditions
Annot.score.up.ancient.mean = aggregate(Module.annot.up.ancient$Jaccard ~ Module.annot.up.ancient$Annot, FUN = mean)
Annot.score.up.ancient.med = aggregate(Module.annot.up.ancient$Jaccard ~ Module.annot.up.ancient$Annot, FUN = median)
Annot.score.up.ancient.sd = aggregate(Module.annot.up.ancient$Jaccard ~ Module.annot.up.ancient$Annot, FUN = sd)
AA = table(Module.annot.up.ancient$Annot)
AA = AA[names(AA) %in% Annot.score.up.ancient.med[,1]]
Annot.score.up.ancient = cbind.data.frame(Annot.score.up.ancient.mean, "Freq" = AA, "median" = Annot.score.up.ancient.med[,2], "SD" = Annot.score.up.ancient.med[,2])
Annot.score.up.ancient = Annot.score.up.ancient[order(Annot.score.up.ancient$median), ]
Annot.score.up.ancient[, 1] = sub("X", "", Annot.score.up.ancient[,1])
Annot.score.up.ancient.3 = Annot.score.up.ancient[Annot.score.up.ancient$Freq > 3,]
Annot.score.up.ancient.3["Rank"] = 1:nrow(Annot.score.up.ancient.3) 

# Ancient anchorpoints with down regulated conditions
Annot.score.down.ancient.mean = aggregate(Module.annot.down.ancient$Jaccard ~ Module.annot.down.ancient$Annot, FUN = mean)
Annot.score.down.ancient.med = aggregate(Module.annot.down.ancient$Jaccard ~ Module.annot.down.ancient$Annot, FUN = median)
Annot.score.down.ancient.sd = aggregate(Module.annot.down.ancient$Jaccard ~ Module.annot.down.ancient$Annot, FUN = sd)
AA = table(Module.annot.down.ancient$Annot)
AA = AA[names(AA) %in% Annot.score.down.ancient.med[,1]]
Annot.score.down.ancient = cbind.data.frame(Annot.score.down.ancient.mean, "Freq" = AA, "median" = Annot.score.down.ancient.med[,2], "SD" = Annot.score.down.ancient.med[,2])
Annot.score.down.ancient = Annot.score.down.ancient[order(Annot.score.down.ancient$median), ]
Annot.score.down.ancient[, 1] = sub("X", "", Annot.score.down.ancient[,1])
Annot.score.down.ancient.3 = Annot.score.down.ancient[Annot.score.down.ancient$Freq > 3,]
Annot.score.down.ancient.3["Rank"] = 1:nrow(Annot.score.down.ancient.3) 
# --------------------------------------------------------------------------------------------------------------------------------------------------


load("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/PO_CORNET_Annotation.RData")


# --------------------------------------------------------------------------------------------------------------------------------------------------
# CORNET annotation text mining
# --------------------------------------------------------------------------------------------------------------------------------------------------
COR.annot = read.delim("/group/biocomp/projects/CORNET/MA_annotation/cornet_allMA_desc_120509.txt")

# Mining PO terms
COR.annot.PO =foreach (j = 1:nrow(COR.annot)) %do% {
	QQ = sapply(unlist(strsplit(as.character(COR.annot[j,6]), "[,]")), strsplit, "[:~]")
	QQ = QQ[sapply(QQ, length) > 0]
	selection = QQ[sapply(QQ, function(x){x[1] == "PO" | x[1] == " PO"})]
	PO = foreach (i = 1:length(selection), .combine = c) %do% {
		unlist(selection[i])[2]
	}
	as.vector(PO)
}

library("memisc")
COR.annot.PO = sapply(COR.annot.PO, trimws)
names(COR.annot.PO) = COR.annot$Name
save(COR.annot.PO, file = "")
save(COR.annot.PO, file = "/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/PO_CORNET_Annotation.RData")
# --------------------------------------------------------------------------------------------------------------------------------------------------
	
Annot.score = list("UP.Recent" = Annot.score.up.recent, "UP.Ancient" = Annot.score.up.ancient, 
	"Down.Recent" = Annot.score.down.recent, "Down.Ancient" = Annot.score.down.ancient)
TOP = function(x){
	conserved = data.frame(table(unlist(COR.annot.PO[names(COR.annot.PO) %in% x[(nrow(x)-20):nrow(x), 1]]))/20)
	conserved = conserved[order(-conserved$Freq), ]
	diverged = data.frame(table(unlist(COR.annot.PO[names(COR.annot.PO) %in% x[1:20, 1]]))/20)
	diverged = diverged[order(-diverged$Freq), ]
	list("Div" = diverged, "Cons" = conserved)
}


CSEA = function(x){
	conserved = data.frame(table(unlist(COR.annot.PO[names(COR.annot.PO) %in% x[(nrow(x)-20):nrow(x), 1]]))/20)
	conserved = conserved[order(-conserved$Freq), ]
	diverged = data.frame(table(unlist(COR.annot.PO[names(COR.annot.PO) %in% x[1:20, 1]]))/20)
	diverged = diverged[order(-diverged$Freq), ]
	list("Div" = diverged, "Cons" = conserved)
}


lapply(Annot.score, TOP)


# --------------------------------------------------------------------------------------------------------------------------------------------------
# 
# --------------------------------------------------------------------------------------------------------------------------------------------------
xvar = list(Annot.score.up.recent.3, Annot.score.up.ancient.3, Annot.score.down.recent.3, Annot.score.down.ancient.3)

Final.Rank = foreach (m = 1:length(xvar)) %do% {
	print(m)
	RAnk = foreach(i = which(sapply(COR.annot.PO[names(COR.annot.PO) %in% xvar[m][[1]][,1]], length) > 0), .combine = rbind.data.frame) %dopar% {
		if(length(unlist(COR.annot.PO[names(COR.annot.PO) %in% xvar[m][[1]][i,1]])) > 0) {
			z = data.frame("PO" = COR.annot.PO[names(COR.annot.PO) %in% xvar[m][[1]][i,1]], "median" = xvar[m][[1]][i,"median"])
			colnames(z) = c("PO", "median")
			z
		}
	}

	Rank.result = as.data.frame(cbind("Median" = aggregate(RAnk$median ~ RAnk$PO, FUN = median)[,2], 
		"Freq" = table(RAnk$PO), "IQR" = aggregate(RAnk$median ~ RAnk$PO, FUN = IQR)[,2], "mean" = aggregate(RAnk$median ~ RAnk$PO, FUN = mean)[,2]))

	Rank.result = Rank.result[order(Rank.result$Median), ]
	Rank.result["PO"] = rownames(Rank.result)

	# --------------------------------------------------------------------------------------------------------------------------------------------------
	# Permutation test
	# --------------------------------------------------------------------------------------------------------------------------------------------------
	PO.size = sort(unique(Rank.result$Freq))
	Permut = matrix(NA, 10000, length(PO.size))
	j = 0
	for (i in PO.size) {
		print(i)
		j = j + 1
		PO.Rand = rep(PO.size[j], 10000)
		# Sample from nodes with the cluster's size
		Permut[, j] = unlist(sapply(PO.Rand, function(x){median(sample(RAnk$median, x))}))
	}

	Thresh = data.frame("Freq" = PO.size, "Threshold.Div" = apply(Permut, 2, quantile, probs = 0.05), 
		"Threshold.Cons" = apply(Permut, 2, quantile, probs = 0.95))
	Rank.result.fin = merge(Thresh, Rank.result, by = "Freq")
	Rank.result.fin.Div = Rank.result.fin[Rank.result.fin$Median < Rank.result.fin$Threshold.Div,]
	Rank.result.fin.Cons = Rank.result.fin[Rank.result.fin$Median > Rank.result.fin$Threshold.Cons,]
	list("Total" = Rank.result.fin, "Diverged" = Rank.result.fin.Div, "Conserved" = Rank.result.fin.Cons)
}
# --------------------------------------------------------------------------------------------------------------------------------------------------
save(Final.Rank, file = "/scratch/iwt/ehsab/Gene_duplication/Results/SA/CORNET_annotaion/Diverged.annot.RData")



