
load("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Compendium_wise/development/Score99/Result/Total1.RData")
A = sapply(g, unlist)
AA = sapply(A, length)
load("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Compendium_wise/development/Score99/Result.RData")
load("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Compendium_wise/development/Score99/G.score.RData")


## Grape
	library(igraph)
	data.grape = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/Vitis_vinifera/GE.mean.SD.txt", nrows = 23629, comment.char = "", header = T)
	gene.names <- read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/Vitis_vinifera/gene.names.txt")
	gene.names = as.matrix(gene.names)
	


	condition.names = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/condition_name.txt", header = T)
	NAME = names(table(condition.names))

	## Arabidopsis
	data.arab = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/Integrated_SD.txt", colClasses = "numeric", nrows = 19285, comment.char = "", header = T)
	data.arab = as.matrix(data.arab[, condition.names == "development"])
	Ref = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/Reference_SD.txt", header = T)
	rownames(data.arab) = Ref[,1]
	
	
	
	k = 6
	
	## Grape modules
	load("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Graph_V2.R")
	g.grape = g
	
	## Arabidopsis modules
	load(paste("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Gene_score/Compendium_wise/",NAME[k],"/Graph_V2.RData", sep = ""))
	
	orthologous = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Orthology/Grape-Arabidopsis/All.txt", header = T)
	onetone = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Orthology/Grape-Arabidopsis/ZhenProcessed.txt", header = T)
	PollenRootLeafclean
	gene.names = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/Vitis_vinifera/gene.names.txt")
	map = read.delim("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Orthology/Grape-Arabidopsis/Map.12x.plaza.txt")
	Ref = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/Reference_SD.txt", header = T)
	# paralogous = read.table("~/Gene_duplication/Seed_genes/paralogous_total.txt", header = T)
	paralogous = read.table("~/Gene_duplication/Seed_genes/paralogous_NoFilter.txt", header = T)
	rownames(data.arab) = Ref[,1]


	

	
library(foreach)
library(doMC)
registerDoMC(40)	

Correlation_total = foreach (m = 1:length(which(AA > 1))) %dopar% {
	print(m)
	LEN = length(unlist(A[which(AA > 1)[m]]))
	Correlation = list()	
	for (n in 1:LEN)	{
		
		i = unlist(A[which(AA > 1)[m]])[n]
		
		memb.comun.grape = data.frame("name" = gene.names[neighborhood(g.grape, 1, i)[[1]],])
		memb.comun.grape.para = merge(memb.comun.grape, map, by = "name")
		## Remove duplicates 
		memb.comun.grape.para = memb.comun.grape.para[!(duplicated(memb.comun.grape.para$name) 
							| duplicated(memb.comun.grape.para$name, fromLast = TRUE)), ]
		memb.comun.grape.para = merge(memb.comun.grape.para, orthologous, by = "PLAZA.ID")
		
		# First pair arabidopsis module to vitis module
		j = Result[which(AA > 1)[m], "Gene"]
		memb.comun.Arab = data.frame("gene_id.y" = Ref[neighborhood(g, 1, j)[[1]],1])
		# Orthologous genes between each module
		Merge1 = merge(memb.comun.Arab, memb.comun.grape.para, by = "gene_id.y")
		Orth1 = Merge1$name
		Orth1.arab = Merge1$gene_id.y
		# Gene expression profile for orthologous genes
		gene.exp1.grape = data.grape[Orth1, ]
		gene.exp1.arab = data.arab[Orth1.arab, ]
		
		# Second pair arabidopsis module to vitis module
		j = Result[which(AA > 1)[m], "Paralogous"]
		memb.comun.Arab = data.frame("gene_id.y" = Ref[neighborhood(g, 1, j)[[1]],1])
		# Orthologous genes between each module
		Merge2 = merge(memb.comun.Arab, memb.comun.grape.para, by = "gene_id.y")
		Orth2 = Merge2$name
		Orth2.arab = Merge2$gene_id.y
		# Gene expression profile for orthologous genes
		gene.exp2.grape = data.grape[Orth2, ]
		gene.exp2.arab = data.arab[Orth2.arab, ]
		
		
		# Whitin correlation between paralogous genes in each module
		cor1.arab = cor(t(gene.exp1.arab))
		cor2.arab = cor(t(gene.exp2.arab))
		cor1.grape = cor(t(gene.exp1.grape))
		cor2.grape = cor(t(gene.exp2.grape))

		COrrelation1 = COrrelation2 = c()
		
		# Between orthologous gene correlation 
		for (i in 1:nrow(cor1.arab)) {
			COrrelation1[i] = cor(cbind(cor1.arab[i,], cor1.grape[i,]))[1,2]
		}
		
		for (i in 1:nrow(cor2.arab)) {
			COrrelation2[i] = cor(cbind(cor2.arab[i,], cor2.grape[i,]))[1,2]
		}
		


		# For mutual orthologous

		Orth = intersect(Orth1, Orth2)
		if (length(Orth) > 1) {

			Orth1.arab.mutual = (Merge1[Merge1$name %in% Orth, "gene_id.y"])
			Orth2.arab.mutual = (Merge2[Merge2$name %in% Orth, "gene_id.y"])

			Orth1.grape.mutual = Merge1[Merge1$gene_id.y %in% Orth1.arab.mutual, "name"]
			Orth2.grape.mutual = Merge2[Merge2$gene_id.y %in% Orth2.arab.mutual, "name"]

			gene.exp1.grape.mutual = data.grape[Orth1.grape.mutual, ]
			gene.exp2.grape.mutual = data.grape[Orth2.grape.mutual, ]
			gene.exp1.arab.mutual = data.arab[Orth1.arab.mutual, ]
			gene.exp2.arab.mutual = data.arab[Orth2.arab.mutual, ]


			cor1.grape.mutual = cor(t(gene.exp1.grape.mutual))
			cor2.grape.mutual = cor(t(gene.exp2.grape.mutual))
			cor1.arab.mutual = cor(t(gene.exp1.arab.mutual))
			cor2.arab.mutual = cor(t(gene.exp2.arab.mutual))
			
			COrrelation2.mutual = COrrelation1.mutual = c()

			for (k in 1:nrow(cor2.arab.mutual)) {
				COrrelation2.mutual[k] = 1/( mean(apply(cbind(cor2.arab.mutual[k,], cor2.grape.mutual[k,]), 1, dist)) + 0.00000001)
			}
			
			for (k in 1:nrow(cor1.arab.mutual)) {
				COrrelation1.mutual[k] = 1/( mean(apply(cbind(cor1.arab.mutual[k,], cor1.grape.mutual[k,]), 1, dist)) + 0.00000001)
			}
		} else {COrrelation1.mutual = COrrelation2.mutual = 0}

		Correlation[n] = list(list(COrrelation1, COrrelation2, COrrelation1.mutual ,COrrelation2.mutual))
		
	}	
	Correlation

}

save(Correlation_total, file = "/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Compendium_wise/development/ICC.RData")

# for (i in 1:LEN){
# 	print(ks.test(unlist(Correlation[i][[1]][1]), unlist(Correlation[i][[1]][2]), 
# 		exact = T)$p)
# }

CORR = foreach (i = 1:length(which(AA > 1)), .combine = rbind) %do% {
	print(i)
	sapply(Correlation_total[[i]][[1]], median)
}
	
	
	
####################################################################################################################################################################### 
#Pollen-Root-Leaf
####################################################################################################################################################################### 

# Between orthologous gene correlation 
CORR = function(x,y) {
	library(foreach)
	AA = foreach (i = 1:nrow(x), .combine = c) %do% {
		cor(cbind(x[i,], y[i,]))[1,2]
	}
	return(AA)
}

PollenRootLeafclean = read.delim("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Orthology/Grape-Arabidopsis/PollenRootLeafcleanGrape12.txt")

Pollen = as.character(PollenRootLeafclean[PollenRootLeafclean$Name %in% jjj & PollenRootLeafclean$Module == "Pollen", "gene_id.y"])
# BPollen = total_gene_scores.a[Root, Root]
BPollen = data.arab[Pollen,]
PollenCor = cor(t(BPollen))


PollenLeaf = as.character(PollenRootLeafclean[PollenRootLeafclean$Name %in% jjj & PollenRootLeafclean$Module == "Pollen", "gene_id.y.1"])
# BPollen = total_gene_scores.a[Root, Root]
BPollenLeaf = data.arab[PollenLeaf,]
PollenLeafCor = cor(t(BPollenLeaf))

BGrape = data.grape[jjj,]

GrapeCor = cor(t(BGrape))
PollenGrape = CORR(GrapeCor, PollenCor)
PollenLeafGrape = CORR(GrapeCor, PollenLeafCor)
RootGrape = CORR(grapeRoot, RootCor)
RootLeafGrape = CORR(grapeRoot, RootLeafCor)


PollenICC = data.frame("Pollen" = Pollen, "Leaf" = PollenLeaf, "PollenICC" = PollenGrape,  "LeafICC" = PollenLeafGrape, "Grape" = jjj)





Root = as.character(PollenRootLeafclean[PollenRootLeafclean$Name %in% jjj & PollenRootLeafclean$Module == "Root", "gene_id.y"])
# BPollen = total_gene_scores.a[Root, Root]
BRoot = data.arab[Root,]
RootCor = cor(t(BRoot))


RootLeaf = as.character(PollenRootLeafclean[PollenRootLeafclean$Name %in% jjj & PollenRootLeafclean$Module == "Root", "gene_id.y.1"])
# BPollen = total_gene_scores.a[Root, Root]
BRootLeaf = data.arab[RootLeaf,]
RootLeafCor = cor(t(BRootLeaf))

BGrape = data.grape[jjj,]

GrapeCor = cor(t(BGrape))
RootGrape = CORR(GrapeCor, RootCor)
RootLeafGrape = CORR(GrapeCor, RootLeafCor)


RootICC = data.frame("Root" = Root, "Leaf" = RootLeaf, "RootICC" = RootGrape,  "LeafICC" = RootLeafGrape, "Grape" = jjj)



		



CORR = function(x,y) {
	library(weights)
	weight = rep(1, nrow(x))
	w = rep(0, nrow(x))
	
	for (i = 1:nrow(x)) {
		weight[i] = cor(cbind(x[i,], y[i,]))[1,2]
		while(abs(w[i] - weight[i]) > 0.001) {
			w[i] = wtd.cors(cbind(x[i,], y[i,]), weight[i])[1,1]
			if(abs(w[i] - weight[i]) > 0.001) 
		}
	}
	return(AA)
}
	
	
##############################################################	
# ICC based on gene scores
##############################################################	

# Arabidopsis
list = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Scripts/General/order.gene.score.txt")
list = as.vector(unlist(list))
setwd("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Gene_score/Compendium_wise/development/Raw/")
G.a = lapply(list, function(x) mget(load(x)))
total_gene_scores.a = matrix(NA, length(unlist(G.a[1])), length(list))

for (i in 1:length(list)){
	print(i)
	total_gene_scores.a[, i] = as.matrix(unlist(G.a[i]))
}

Ref = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/Reference_SD.txt", header = T)
colnames(total_gene_scores.a) = rownames(total_gene_scores.a) = Ref[, 1]

# Grape
list = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Scripts/General/order.gene_score.Grape.txt")
list = as.vector(unlist(list))
setwd("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Raw/Gene_score/")
G = lapply(list, read.table, header = T)
total_gene_scores = matrix(NA, length(unlist(G[1])), length(list))

for (i in 1:length(list)){
	print(i)
	total_gene_scores[, i] = as.matrix(unlist(G[i]))
}

colnames(total_gene_scores) = rownames(total_gene_scores) = gene.names[, 1]

	
##################################################################################################################################################	
# Comparison of Gene Coexpression Profiles and Construction of Conserved Gene Networks to Find Functional Modules
# Kinoshita
##################################################################################################################################################		
library(foreach)
library(doMC)
registerDoMC(5)
options(width = 200)
library(MESS)

gene.names <- read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/Vitis_vinifera/gene.names.txt")
gene.names = as.matrix(gene.names)

orthologous = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Orthology/Grape-Arabidopsis/All.txt", header = T)

###############################################################################################################################################################
# Remove non-uniqued mapped PLAZA ID to gene names and the other way around
###############################################################################################################################################################
map = read.delim("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Orthology/Grape-Arabidopsis/Map.12x.plaza.txt")
# selection = which(sapply(1:sum(duplicated(map$name)), function(i){length(unique(map[map$name == map$name[duplicated(map$name)][i],"PLAZA.ID"]))}) > 1)
# map = map[!(map$name %in% map$name[duplicated(map$name)][selection]),]
# selection = which(sapply(1:sum(duplicated(map$PLAZA.ID)), function(i){length(unique(map[map$PLAZA.ID == map$PLAZA.ID[duplicated(map$PLAZA.ID)][i],"name"]))}) > 1)
# map = map[!(map$PLAZA.ID %in% map$PLAZA.ID[duplicated(map$PLAZA.ID)][selection]),]
map = map[!(duplicated(map[, c("name", "PLAZA.ID")])),]

###############################################################################################################################################################

Ref = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/Reference_SD.txt", header = T)
rownames(Ref) = Ref[,1]
onetone = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Orthology/Grape-Arabidopsis/ZhenProcessed.txt", header = T)
onetone = onetone[!(duplicated(onetone[, "PLAZA.ID"]) | duplicated(onetone[, "gene_id.y"])),]
map1 = merge(map, onetone, by = "PLAZA.ID")


		
OrthZhen = read.delim("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Orthology/Grape-Arabidopsis/PollenRootLeafclean.txt")
Pollen = as.character(OrthZhen[OrthZhen$Module == "Pollen", "gene_id.y",])
PollenLeaf = as.character(OrthZhen[OrthZhen$Module == "Pollen", "gene_id.y.1",])
PollenGrape = as.character(OrthZhen[OrthZhen$Module == "Pollen", "name",])

Root = as.character(OrthZhen[OrthZhen$Module == "Root", "gene_id.y",])
RootLeaf = as.character(OrthZhen[OrthZhen$Module == "Root", "gene_id.y.1",])
RootGrape = as.character(OrthZhen[OrthZhen$Module == "Root", "name",])


PolleN = foreach (k = 1:length(Pollen), .combine = rbind) %dopar% {
	print(k)
	load(paste("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Gene_score/Compendium_wise/development/Raw/G",Ref[Pollen[k], 2],".RData", sep = ""))
	gene_scorespollen = gene_scores
	rownames(gene_scorespollen) = Ref[,1]
	gene_scorespollen = gene_scorespollen[rownames(gene_scorespollen) %in% unique(as.character(map1$gene_id.y)),]
	gene_scorespollenUp = gene_scorespollen[order(-gene_scorespollen)]
	gene_scorespollenDown = gene_scorespollen[order(gene_scorespollen)]

	load(paste("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Gene_score/Compendium_wise/development/Raw/G",Ref[PollenLeaf[k], 2],".RData", sep = ""))
	gene_scorespollenleaf = gene_scores
	rownames(gene_scorespollenleaf) = Ref[,1]
	gene_scorespollenleaf = gene_scorespollenleaf[rownames(gene_scorespollenleaf) %in% unique(as.character(map1$gene_id.y)),]
	gene_scorespollenleafUp = gene_scorespollenleaf[order(-gene_scorespollenleaf)]
	gene_scorespollenleafDown = gene_scorespollenleaf[order(gene_scorespollenleaf)]

	gene_scoregrape = as.matrix(read.delim(paste("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Raw/Gene_score/G_Total.SD",which(gene.names %in% PollenGrape[k]),".txt", sep = "")))
	rownames(gene_scoregrape) = gene.names[,1]
	gene_scoregrape = gene_scoregrape[rownames(gene_scoregrape) %in% unique(as.character(map1$name)),]
	gene_scoregrapeUp = gene_scoregrape[order(-gene_scoregrape)]
	gene_scoregrapeDown = gene_scoregrape[order(gene_scoregrape)]

	Ranks = merge(merge(merge(
	map1,
	data.frame("gene_id.y" = names(gene_scorespollenUp), "RankPollen" = 1:length(gene_scorespollenUp))),
	data.frame("gene_id.y" = names(gene_scorespollenleafUp), "RankLeaf" = 1:length(gene_scorespollenleafUp))),
	data.frame("name" = names(gene_scoregrapeUp), "RankGrape" = 1:length(gene_scoregrapeUp)))
	Ranks$RankPollen = rank(Ranks$RankPollen)
	Ranks$RankLeaf = rank(Ranks$RankLeaf)
	Ranks$RankGrape = rank(Ranks$RankGrape)
	Ronetone = data.frame("AllRanksPollen" = median(abs(Ranks$RankPollen - Ranks$RankGrape)), "AllRanksLeaf" = median(abs(Ranks$RankLeaf - Ranks$RankGrape)), "Pvalue" = wilcox.test(abs(Ranks$RankPollen - Ranks$RankGrape), abs(Ranks$RankLeaf - Ranks$RankGrape))$p.value)


	AA = foreach (i = 1:500, .combine = rbind) %do% {
		# print(i)
		data.frame("N" = i,	"PollenScore" = length(intersect(as.character(unique(map1[map1$gene_id.y %in% names(gene_scorespollenUp)[1:i], "name"])), names(gene_scoregrapeUp)[1:i])),
			"LeafScore" = length(intersect(as.character(unique(map1[map1$gene_id.y %in% names(gene_scorespollenleafUp)[1:i], "name"])), names(gene_scoregrapeUp)[1:i])))
	}

	BB = foreach (i = 1:500, .combine = rbind) %do% {
		# print(i)
		data.frame("N" = i,	"PollenScore" = length(intersect(as.character(unique(map1[map1$gene_id.y %in% names(gene_scorespollenDown)[1:i], "name"])), names(gene_scoregrapeDown)[1:i])),
			"LeafScore" = length(intersect(as.character(unique(map1[map1$gene_id.y %in% names(gene_scorespollenleafDown)[1:i], "name"])), names(gene_scoregrapeDown)[1:i])))
	}

	DD = cbind("N" = AA$N, AA[,2:3] + BB[,2:3])

	data.frame("Pollen" = Pollen[k], "Leaf" = PollenLeaf[k], Ronetone,  "Grape" = PollenGrape[k],  DD[nrow(DD),2:3], PollenAUC = auc(DD$N, DD$PollenScore)/auc(1:nrow(DD), 1:nrow(DD))
		, LeafAUC = auc(DD$N, DD$LeafScore)/auc(1:nrow(DD), 1:nrow(DD)), "UpPollen" = AA[nrow(AA),2], "UpLeaf" = AA[nrow(AA),3], "DownPollen" = BB[nrow(BB),2], "DownLeaf" = BB[nrow(BB),3])
}



RooT = foreach (k = 1:length(Root), .combine = rbind) %dopar% {
	print(k)
	load(paste("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Gene_score/Compendium_wise/development/Raw/G",Ref[Root[k], 2],".RData", sep = ""))
	gene_scoresRoot = gene_scores
	rownames(gene_scoresRoot) = Ref[,1]
	gene_scoresRoot = gene_scoresRoot[rownames(gene_scoresRoot) %in% unique(as.character(map1$gene_id.y)),]
	gene_scoresRootUp = gene_scoresRoot[order(-gene_scoresRoot)]
	gene_scoresRootDown = gene_scoresRoot[order(gene_scoresRoot)]

	load(paste("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Gene_score/Compendium_wise/development/Raw/G",Ref[RootLeaf[k], 2],".RData", sep = ""))
	gene_scoresRootleaf = gene_scores
	rownames(gene_scoresRootleaf) = Ref[,1]
	gene_scoresRootleaf = gene_scoresRootleaf[rownames(gene_scoresRootleaf) %in% unique(as.character(map1$gene_id.y)),]
	gene_scoresRootleafUp = gene_scoresRootleaf[order(-gene_scoresRootleaf)]
	gene_scoresRootleafDown = gene_scoresRootleaf[order(gene_scoresRootleaf)]

	gene_scoregrape = as.matrix(read.delim(paste("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Raw/Gene_score/G_Total.SD",which(gene.names %in% RootGrape[k]),".txt", sep = "")))
	rownames(gene_scoregrape) = gene.names[,1]
	gene_scoregrape = gene_scoregrape[rownames(gene_scoregrape) %in% unique(as.character(map1$name)),]
	gene_scoregrapeUp = gene_scoregrape[order(-gene_scoregrape)]
	gene_scoregrapeDown = gene_scoregrape[order(gene_scoregrape)]


	Ranks = merge(merge(merge(
	map1,
	data.frame("gene_id.y" = names(gene_scoresRootUp), "RankRoot" = 1:length(gene_scoresRootUp))),
	data.frame("gene_id.y" = names(gene_scoresRootleafUp), "RankLeaf" = 1:length(gene_scoresRootleafUp))),
	data.frame("name" = names(gene_scoregrapeUp), "RankGrape" = 1:length(gene_scoregrapeUp)))
	Ranks$RankRoot = rank(Ranks$RankRoot)
	Ranks$RankLeaf = rank(Ranks$RankLeaf)
	Ranks$RankGrape = rank(Ranks$RankGrape)
	Ronetone = data.frame("AllRanksRoot" = median(abs(Ranks$RankRoot - Ranks$RankGrape)), "AllRanksLeaf" = median(abs(Ranks$RankLeaf - Ranks$RankGrape)), "Pvalue" = wilcox.test(abs(Ranks$RankRoot - Ranks$RankGrape), abs(Ranks$RankLeaf - Ranks$RankGrape))$p.value)


	AA = foreach (i = 1:500, .combine = rbind) %do% {
		# print(i)
		data.frame("N" = i,	"RootScore" = length(intersect(as.character(unique(map1[map1$gene_id.y %in% names(gene_scoresRootUp)[1:i], "name"])), names(gene_scoregrapeUp)[1:i])),
			"LeafScore" = length(intersect(as.character(unique(map1[map1$gene_id.y %in% names(gene_scoresRootleafUp)[1:i], "name"])), names(gene_scoregrapeUp)[1:i])))
	}

	BB = foreach (i = 1:500, .combine = rbind) %do% {
		# print(i)
		data.frame("N" = i,	"RootScore" = length(intersect(as.character(unique(map1[map1$gene_id.y %in% names(gene_scoresRootDown)[1:i], "name"])), names(gene_scoregrapeDown)[1:i])),
			"LeafScore" = length(intersect(as.character(unique(map1[map1$gene_id.y %in% names(gene_scoresRootleafDown)[1:i], "name"])), names(gene_scoregrapeDown)[1:i])))
	}

	DD = cbind("N" = AA$N, AA[,2:3] + BB[,2:3])

	data.frame("Root" = Root[k], "Leaf" = RootLeaf[k], Ronetone, "Grape" = RootGrape[k], DD[nrow(DD),2:3], RootAUC = auc(DD$N, DD$RootScore)/auc(1:nrow(DD), 1:nrow(DD))
		, LeafAUC = auc(DD$N, DD$LeafScore)/auc(1:nrow(DD), 1:nrow(DD)), "UpRoot" = AA[nrow(AA),2], "UpLeaf" = AA[nrow(AA),3], "DownRoot" = BB[nrow(BB),2], "DownLeaf" = BB[nrow(BB),3])
}


PolleN$ScoreFDR = p.adjust(sapply(1:nrow(PolleN), function(i){prop.test(x = c(PolleN$PollenScore[i], PolleN$LeafScore[i]), n = c(1000, 1000))$p.value}), method = "fdr")
RooT$ScoreFDR = p.adjust(sapply(1:nrow(RooT), function(i){prop.test(x = c(RooT$RootScore[i], RooT$LeafScore[i]), n = c(1000, 1000))$p.value}), method = "fdr")


# RooT$Module = ifelse(RooT$AllRanksRoot < RooT$AllRanksLeaf & RooT$Pvalue < 0.05, "Root", "Leaf")
RooT$Module = ifelse(RooT$RootScore > RooT$LeafScore & RooT$ScoreFDR < 0.1, "Root", ifelse(RooT$RootScore < RooT$LeafScore & RooT$ScoreFDR < 0.1, "Leaf", "Both"))
# PolleN$Module = ifelse(PolleN$AllRanksPollen < PolleN$AllRanksLeaf & PolleN$Pvalue < 0.05, "Pollen", "Leaf")
PolleN$Module = ifelse(PolleN$PollenScore > PolleN$LeafScore & PolleN$ScoreFDR < 0.1, "Pollen", ifelse(PolleN$PollenScore < PolleN$LeafScore & PolleN$ScoreFDR < 0.1, "Leaf", "Both"))



write.table(PolleN, file = "/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Grape/Multiple_othologous_modules/Pollen.txt", sep = "\t", row.names = F, quote = F)
write.table(RooT, file = "/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Grape/Multiple_othologous_modules/Root.txt", sep = "\t", row.names = F, quote = F)



##################################################################################################################################################	
# ICC for one to one orthologes in each module with the Grape 
##################################################################################################################################################		

# Correlation between two rows of two different matrices
source("/ngsprojects/iwt/ehsab/Gene_duplication/Scripts/Functions/CORR.R")
Ref = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/Reference_SD.txt", header = T)
rownames(Ref) = Ref[,1]


LeafRoot = read.table("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/Module/AT3G19710_Leaf_90.txt")
Root = read.table("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/Module/AT1G50090_Root_90.txt")
AnchRoot = read.delim("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/Anchorpoint/AT1G50090_AT3G19710_52_90.txt")



Leaf = read.table("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/Module/Leaf_90.txt")
Pollen = read.table("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/Module/Pollen_90.txt")
Anch = read.delim("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/Anchorpoint/Geneset1_2_53_90.txt.txt")

map = read.delim("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Orthology/Grape-Arabidopsis/Map.12x.plaza.txt")
map = map[!(duplicated(map[, c("name", "PLAZA.ID")])),] # Map PLAZA ID to the gene names in the grapevine expression data
onetone = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Orthology/Grape-Arabidopsis/ZhenProcessed.txt", header = T)
maponetone = merge(map, onetone, by = "PLAZA.ID")
maponetone = maponetone[maponetone$gene_id.y %in% Ref$Gene,] # Remove arabidopsis genes which are not in the compendium
maponetone = maponetone[maponetone$name %in% rownames(data.grape),] # Remoe grapevine genes which are not in the development compendium




require(caret)
flds = c(createFolds(1:ncol(data.grape), k = 10, list = TRUE, returnTrain = FALSE),
	 createFolds(1:ncol(data.grape), k = 10, list = TRUE, returnTrain = FALSE),
	 createFolds(1:ncol(data.grape), k = 10, list = TRUE, returnTrain = FALSE),
	 createFolds(1:ncol(data.grape), k = 10, list = TRUE, returnTrain = FALSE),
	 createFolds(1:ncol(data.grape), k = 10, list = TRUE, returnTrain = FALSE),
	 createFolds(1:ncol(data.grape), k = 10, list = TRUE, returnTrain = FALSE),
	 createFolds(1:ncol(data.grape), k = 10, list = TRUE, returnTrain = FALSE),
	 createFolds(1:ncol(data.grape), k = 10, list = TRUE, returnTrain = FALSE),
	 createFolds(1:ncol(data.grape), k = 10, list = TRUE, returnTrain = FALSE),
	 createFolds(1:ncol(data.grape), k = 10, list = TRUE, returnTrain = FALSE))

RR = foreach(k = 1:length(flds)) %dopar% {
	print(k)
	COR.grape = cor(t(data.grape[as.character(maponetone$name), flds[k][[1]]]))
	COR.arab = cor(t(data.arab[as.character(maponetone$gene_id.y),]))
	onetoneCorr = CORR(COR.grape, COR.arab, w = runif(nrow(COR.grape)))

	threshold = 1
	i = 0

	while (threshold > 0.1) {
		i = i +1
		# print(i)
		Exclude = which(onetoneCorr < 0)
		COR.grape = COR.grape[-Exclude,-Exclude]
		COR.arab = COR.arab[-Exclude,-Exclude]
		onetoneCorrold = onetoneCorr[-Exclude]
		onetoneCorr = CORR(COR.grape, COR.arab, w = onetoneCorrold)
		threshold = sum((onetoneCorr - onetoneCorrold)^2)
	}

	save(COR.grape, COR.arab, onetoneCorr, file = paste("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Permutation/10Fold/",k,".RData", sep = ""))

}



	COR.grape = cor(t(data.grape[as.character(maponetone$name), ]))
	COR.arab = cor(t(data.arab[as.character(maponetone$gene_id.y),]))
	onetoneCorr = CORR(COR.grape, COR.arab, w = rep(1, nrow(COR.grape)))

	threshold = 1
	i = 0

	while (threshold > 0.1) {
		i = i +1
		# print(i)
		Exclude = which(onetoneCorr < 0)
		COR.grape = COR.grape[-Exclude,-Exclude]
		COR.arab = COR.arab[-Exclude,-Exclude]
		onetoneCorrold = onetoneCorr[-Exclude]
		onetoneCorr = CORR(COR.grape, COR.arab, w = onetoneCorrold)
		threshold = sum((onetoneCorr - onetoneCorrold)^2)
	}




onetoneCorr1 = onetoneCorr

COR.arab = RR[i][[1]][[2]]
onetoneCorr = RR[i][[1]][[3]]

lapply(1:10, function(i){quantile(RR[i][[1]][[3]][which(rownames(RR[i][[1]][[2]]) %in% setdiff(LeafRoot[,1],intersect(LeafRoot[,1], Leaf[,1])))])})
lapply(1:10, function(i){quantile(RR[i][[1]][[3]][which(rownames(RR[i][[1]][[2]]) %in% setdiff(Leaf[,1],intersect(LeafRoot[,1], Leaf[,1])))])})
lapply(1:10, function(i){quantile(RR[i][[1]][[3]][which(rownames(RR[i][[1]][[2]]) %in% intersect(Leaf[,1], LeafRoot[,1]))])})
lapply(1:10, function(i){quantile(RR[i][[1]][[3]][which(rownames(RR[i][[1]][[2]]) %in% setdiff(Root[,1],intersect(Pollen[,1], Root[,1])))])})
lapply(1:10, function(i){quantile(RR[i][[1]][[3]][which(rownames(RR[i][[1]][[2]]) %in% setdiff(Pollen[,1],intersect(Pollen[,1], Root[,1])))])})
lapply(1:10, function(i){quantile(RR[i][[1]][[3]][which(rownames(RR[i][[1]][[2]]) %in% intersect(Pollen[,1], Root[,1]))])})
lapply(1:10, function(i){quantile(RR[i][[1]][[3]][which(rownames(RR[i][[1]][[2]]) %in% Pollen[,1])])})
lapply(1:10, function(i){quantile(RR[i][[1]][[3]][which(rownames(RR[i][[1]][[2]]) %in% Root[,1])])})
lapply(1:10, function(i){quantile(RR[i][[1]][[3]][which(rownames(RR[i][[1]][[2]]) %in% LeafRoot[,1])])})



RR = foreach(k = 1:100, .combine = rbind) %dopar% {
	print(k)
	load(paste("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Permutation/10Fold/",k,".RData", sep = ""))
	cbind.data.frame(rbind(cbind.data.frame(t(quantile(onetoneCorr[which(rownames(COR.arab) %in% Root[,1])])), "Module" = "Root"),
	cbind.data.frame(t(quantile(onetoneCorr[which(rownames(COR.arab) %in% LeafRoot[,1])])), "Module" = "LeafRoot"),
	cbind.data.frame(t(quantile(onetoneCorr[which(rownames(COR.arab) %in% Pollen[,1])])), "Module" = "Pollen"),
	cbind.data.frame(t(quantile(onetoneCorr[which(rownames(COR.arab) %in% Leaf[,1])])), "Module" = "LeafPollen"),
	cbind.data.frame(t(quantile(onetoneCorr[which(rownames(COR.arab) %in% setdiff(LeafRoot[,1],intersect(LeafRoot[,1], Leaf[,1])))])), "Module" = "LeafRootS"),
	cbind.data.frame(t(quantile(onetoneCorr[which(rownames(COR.arab) %in% setdiff(Leaf[,1],intersect(LeafRoot[,1], Leaf[,1])))])), "Module" = "LeafPollenS"),
	cbind.data.frame(t(quantile(onetoneCorr[which(rownames(COR.arab) %in% intersect(Leaf[,1], LeafRoot[,1]))])), "Module" = "LeafLeaf"),
	cbind.data.frame(t(quantile(onetoneCorr[which(rownames(COR.arab) %in% setdiff(Root[,1],intersect(Pollen[,1], Root[,1])))])), "Module" = "RootS"),
	cbind.data.frame(t(quantile(onetoneCorr[which(rownames(COR.arab) %in% setdiff(Pollen[,1],intersect(Pollen[,1], Root[,1])))])), "Module" = "PollenS"),
	cbind.data.frame(t(quantile(onetoneCorr[which(rownames(COR.arab) %in% intersect(Pollen[,1], Root[,1]))])), "Module" = "RootPollen")), "Fold" = k)
}





Graph = rbind(data.frame("x" = onetoneCorr[which(rownames(COR.arab) %in% Root[,1])], "Module" = "Root"),
data.frame("x" = onetoneCorr[which(rownames(COR.arab) %in% LeafRoot[,1])], "Module" = "Leaf-Root"),
data.frame("x" = onetoneCorr[which(rownames(COR.arab) %in% Pollen[,1])], "Module" = "Pollen"),
data.frame("x" = onetoneCorr[which(rownames(COR.arab) %in% Leaf[,1])], "Module" = "Leaf-Pollen")) 


ggsave(file = , plot = ggplot(Graph, aes(x, fill = Module)) + geom_histogram(alpha = 1, binwidth = 0.03))



maponetonePollen = as.character(maponetone[as.character(maponetone$gene_id.y) %in% Pollen[,1], "gene_id.y"])
maponetonePollenGrape = as.character(maponetone[as.character(maponetone$gene_id.y) %in% Pollen[,1], "name"])
maponetonePollenLeaf = as.character(maponetone[as.character(maponetone$gene_id.y) %in% Leaf[,1], "gene_id.y"])
maponetonePollenLeafGrape = as.character(maponetone[as.character(maponetone$gene_id.y) %in% Leaf[,1], "name"])

maponetoneRoot = as.character(maponetone[as.character(maponetone$gene_id.y) %in% Root[,1], "gene_id.y"])
maponetoneRootGrape = as.character(maponetone[as.character(maponetone$gene_id.y) %in% Root[,1], "name"])
maponetoneRootLeaf = as.character(maponetone[as.character(maponetone$gene_id.y) %in% LeafRoot[,1], "gene_id.y"])
maponetoneRootLeafGrape = as.character(maponetone[as.character(maponetone$gene_id.y) %in% LeafRoot[,1], "name"])

maponetoneAll = rbind(data.frame("ORTH" = maponetonePollen, "Module" = "Pollen"),
	data.frame("ORTH" = maponetonePollenGrape, "Module" = "PollenGrape"),
	data.frame("ORTH" = maponetonePollenLeaf, "Module" = "PollenLeaf"),
	data.frame("ORTH" = maponetonePollenLeafGrape, "Module" = "PollenLeafGrape"),
	data.frame("ORTH" = maponetoneRoot, "Module" = "Root"),
	data.frame("ORTH" = maponetoneRootGrape, "Module" = "RootGrape"),
	data.frame("ORTH" = maponetoneRootLeaf, "Module" = "RootLeaf"),
	data.frame("ORTH" = maponetoneRootLeafGrape, "Module" = "RootLeafGrape"))

write.table(maponetoneAll, row.names = F, quote = F, sep = "\t", file = "/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Orthology/Grape-Arabidopsis/ZhenOneTOneModule.txt")

COR.grape.Pollen = cor(t(data.grape[as.character(maponetonePollenGrape),]))
COR.arab.Pollen = cor(t(data.arab[as.character(maponetonePollen),]))
onetoneCorrPollen = CORR(COR.grape.Pollen, COR.arab.Pollen)
COR.grape.PollenLeaf = cor(t(data.grape[as.character(maponetonePollenLeafGrape),]))
COR.arab.PollenLeaf = cor(t(data.arab[as.character(maponetonePollenLeaf),]))
onetoneCorrPollenLeaf = CORR(COR.grape.PollenLeaf, COR.arab.PollenLeaf)


COR.grape.Root = cor(t(data.grape[as.character(maponetoneRootGrape),]))
COR.arab.Root = cor(t(data.arab[as.character(maponetoneRoot),]))
onetoneCorrRoot = CORR(COR.grape.Root, COR.arab.Root)
COR.grape.RootLeaf = cor(t(data.grape[as.character(maponetoneRootLeafGrape),]))
COR.arab.RootLeaf = cor(t(data.arab[as.character(maponetoneRootLeaf),]))
onetoneCorrRootLeaf = CORR(COR.grape.RootLeaf, COR.arab.RootLeaf)

RootPlot = data.frame("ICC" = c(onetoneCorrRoot, onetoneCorrRootLeaf), "Module" = c(rep("Root", length(onetoneCorrRoot)), rep("Leaf", length(onetoneCorrRootLeaf))))
PollenPlot = data.frame("ICC" = c(onetoneCorrPollen, onetoneCorrPollenLeaf), "Module" = c(rep("Pollen", length(onetoneCorrPollen)), rep("Leaf", length(onetoneCorrPollenLeaf))))

ggplot(RootPlot, aes(ICC, fill = Module)) + geom_density(alpha = 0.3)
ggplot(PollenPlot, aes(ICC, fill = Module)) + geom_density(alpha = 0.3)
