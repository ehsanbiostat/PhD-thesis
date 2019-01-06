library(igraph)
library(foreach)
library(doMC)
# registerDoMC(2)
k = 6
condition.names = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated CORNET2.0/condition_name.txt", header = T)
NAME = names(table(condition.names))
NAME.spe = NAME[k]
source("/ngsprojects/iwt/ehsab/Gene_duplication/Scripts/General/myImagePlot.R")
source("/ngsprojects/iwt/ehsab/Gene_duplication/Scripts/General/myImagePlot_sym.R")
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

data.exp = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated CORNET2.0/Integrated_SD.txt", colClasses = "numeric", nrows = 19285, comment.char = "", header = T)
data.exp.scale = apply(data.exp,1,scale)
data.exp.scale = t(data.exp.scale)
data.exp.scale = apply(data.exp.scale, 2, scale)
data.exp.scale = data.exp.scale[, condition.names == NAME.spe]


data.exp.scale.med = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated CORNET2.0/Integrated_SD_Norm_Med.txt", colClasses = "numeric", nrows = 19285, comment.char = "", header = T)

data.exp.scale.med = data.exp.scale.med[, condition.names == NAME.spe]

Sig = ifelse(data.exp.scale > 0, 1, ifelse(data.exp.scale < 0, -1, 0))
Sig.med = ifelse(data.exp.scale.med > 0, 1, ifelse(data.exp.scale.med < 0, -1, 0))

Sig.med.zero = ifelse(data.exp.scale.med > 0, 1, ifelse(data.exp.scale.med < 0, 0, 0))


ref = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated CORNET2.0/Reference_SD.txt", header = T)
paralogous = read.table("~/Gene_duplication/Seed_genes/paralogous_total.txt", header = T)
## Load generated graph from SA result
load(paste("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Gene_score/Compendium_wise/",NAME.spe,"/Graph_V2.RData", sep = ""))


ORTHO = read.table("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Validation/Compendium_wise/development/5.Para.list.txt", header = T)


foreach (N = 1:200) %dopar% {
	print(N)
	j = ORTHO$orthologous[N]
	jpar = ORTHO$orthologous.1[N]
	memb.comun = neighborhood(g,1,j)[[1]]
	memb.comun.para = neighborhood(g,1,jpar)[[1]]
	CondRank1 = foreach(i = 1:length(memb.comun), .combine = rbind) %do%{
		rank(total_cond_scores[, memb.comun[i]])
	}
	CondRank2 = foreach(i = 1:length(memb.comun.para), .combine = rbind) %do%{
		rank(total_cond_scores[, memb.comun.para[i]])
	}
	rownames(CondRank2) = c(1:length(memb.comun.para))
	rownames(CondRank1) = c(1:length(memb.comun))

	Med.cond.score1 = apply(CondRank1, 2, median)
	Med.cond.score2 = apply(CondRank2, 2, median)

	CondRank = rbind(c(1:dim(total_cond_scores)[1]), CondRank1, CondRank2)
	colnames(CondRank) = colnames(data.exp.scale)
	CondRank0 = CondRank[,order(Med.cond.score1, decreasing = T)]
	CondRank = rbind(c(1:dim(total_cond_scores)[1]), CondRank2, CondRank1)
	colnames(CondRank) = colnames(data.exp.scale)
	CondRank01 = CondRank[,order(Med.cond.score2, decreasing = T)]
	# myImagePlot_sym(data.exp.scale[c(memb.comun, memb.comun.para),CondRank0[1,]])
	# myImagePlot_sym(data.exp.scale[c(memb.comun, memb.comun.para),CondRank01[1,]])
	# myImagePlot_sym(data.exp.scale[c(memb.comun, memb.comun.para),c(CondRank0[1,1:20],CondRank01[1,1:20])])

	Inter1 = intersect(memb.comun.para, paralogous[paralogous$Ref %in% memb.comun, "Anch"])
	Inter2 = intersect(memb.comun, paralogous[paralogous$Ref %in% memb.comun.para, "Anch"])
	aa = paralogous[paralogous$Ref %in% Inter2, ]
	bb = aa[aa$Anch %in% Inter1,]
	colnames(bb)[1] = "References"
	Inter_name = merge(ref, bb, by = "References")
	colnames(Inter_name)[c(1,3)] = c("Ref", "References")
	Inter_name = merge(ref, Inter_name, by = "References")
	i = paste(j, jpar, sep = "_")
	write.table(Inter_name, file = paste("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Local/Compendium_wise/",NAME.spe,"/Score99/Module_Anchorpoints_V2/Module",i,".txt", sep = ""), row.names = F, quote = F, sep = "\t")
	k = i
	Result = read.table(paste("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Local/Compendium_wise/",NAME.spe,"/Score99/Module_Anchorpoints_V2/Module",k,".txt", sep = ""), header = T)
	j = Result$Ref #Gene.y
	jpar = Result$References #Gene.x
	
	
	CondRank1 = foreach(i = 1:length(j), .combine = rbind) %do%{
		rank(total_cond_scores[, j[i]])
	}
	CondRank2 = foreach(i = 1:length(jpar), .combine = rbind) %do%{
		rank(total_cond_scores[, jpar[i]])
	}
	
	Med.cond.score1 = apply(CondRank1, 2, median)
	Med.cond.score2 = apply(CondRank2, 2, median)

	CondRank = rbind(c(1:dim(total_cond_scores)[1]), Med.cond.score1, CondRank1, CondRank2)
	colnames(CondRank) = colnames(data.exp.scale.med)
	CondRank0 = CondRank[,order(CondRank["Med.cond.score1", ], decreasing = T)]

	CondRank = rbind(c(1:dim(total_cond_scores)[1]), Med.cond.score2, CondRank2, CondRank1)
	colnames(CondRank) = colnames(data.exp.scale.med)
	CondRank01 = CondRank[,order(CondRank["Med.cond.score2", ], decreasing = T)]

	jpeg(paste("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Validation/Compendium_wise/development/Heatmap/",k,"_Sig_4.jpeg", sep = ""), width = 11, height = 7, units = "in", res = 200, pointsize=8)
	myImagePlot_sym(Sig.med[c(j, jpar),c(CondRank0[1,1:15],CondRank01[1,1:15])])
	dev.off()



A = Sig.med.zero[c(j),c(CondRank0[1,1:10])]
A.par = Sig.med.zero[c(jpar),c(CondRank0[1,1:10])]
Ascore = sum(A)/(nrow(A) * ncol(A))
Ascore.par = sum(A.par)/(nrow(A) * ncol(A))

B = Sig.med.zero[c(j), c(CondRank01[1,1:10])]
B.par = Sig.med.zero[c(jpar), c(CondRank01[1,1:10])]
Bscore = sum(B)/(nrow(B) * ncol(B))
Bscore.par = sum(B.par)/(nrow(B) * ncol(B))

c((Ascore - Ascore.par), (Bscore.par - Bscore))

}

myImagePlot_sym(data.exp.scale[c(j, jpar),CondRank0[1,]])
myImagePlot_sym(data.exp.scale[c(j, jpar),CondRank0[1,1:20]])
myImagePlot_sym(data.exp.scale[c(j, jpar),CondRank01[1,1:20]])


# Make heatmap for triangile genes in the modules
load("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Validation/Compendium_wise/development/5.Para.Orth.Heatmap.R")
foreach (N = 1:length(Inter)) %dopar% {
	print(N)
	j = as.data.frame(Inter[N])[, 1] #Gene.y
	jpar = as.data.frame(Inter[N])[, 2] #Gene.x
	
	j = ref[ref[,1] %in% j, "References"]
	jpar = ref[ref[,1] %in% jpar, "References"]
	

	CondRank1 = foreach(i = 1:length(j), .combine = rbind) %do%{
		rank(total_cond_scores[, j[i]])
	}
	CondRank2 = foreach(i = 1:length(jpar), .combine = rbind) %do%{
		rank(total_cond_scores[, jpar[i]])
	}
	
	Med.cond.score1 = apply(CondRank1, 2, median)
	Med.cond.score2 = apply(CondRank2, 2, median)

	CondRank = rbind(c(1:dim(total_cond_scores)[1]), Med.cond.score1, CondRank1, CondRank2)
	colnames(CondRank) = colnames(data.exp)
	CondRank0 = CondRank[,order(CondRank["Med.cond.score1", ], decreasing = T)]

	CondRank = rbind(c(1:dim(total_cond_scores)[1]), Med.cond.score2, CondRank2, CondRank1)
	colnames(CondRank) = colnames(data.exp)
	CondRank01 = CondRank[,order(CondRank["Med.cond.score2", ], decreasing = T)]



# ---------------------------------------------------------------------------------------
# Grape
# development
data.grape = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated CORNET2.0/Vitis_vinifera/GE.mean.SD.txt", nrows = 19589, comment.char = "", header = T)
gene.names <- read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated CORNET2.0/Vitis_vinifera/gene.names.txt")
gene.names = as.matrix(gene.names)

# Stress
# data.grape = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated CORNET2.0/Vitis_vinifera/Stress/Exp.matrix.txt", nrows = 7179, comment.char = "", header = T)
# gene.names <- rownames(data.grape)
gene.names = as.matrix(gene.names)
data.grape.scale = apply(data.grape, 1, scale)
data.grape.scale = t(data.grape.scale)

# development
load("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Graph_V2.R")

# Stress
# load("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Compendium_wise/stress/Graph_V2.R")
	g.grape = g
	
list = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Scripts/General/order.condition_score.Grape.txt")
list = as.vector(unlist(list))
## Read text or RData file with read.table or load function
# development
setwd("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Raw/Condition_score/")
# Stress
# setwd("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Compendium_wise/stress/Raw/Condition_score/")

# C = lapply(list, read.table, sep= "\t" , header=T)
C = lapply(list, read.table, header = T)
total_cond_scores = matrix(NA, length(unlist(C[1])), length(list))
for (i in 1:length(list)){
print(i)
total_cond_scores[, i] = as.matrix(unlist(C[i]))
}
colnames(total_cond_scores) = gene.names[, 1]



orthologous = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Orthology/Grape-Arabidopsis/All.txt", header = T)
map = read.delim("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Orthology/Grape-Arabidopsis/Map.12x.plaza.txt")
Ref = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated CORNET2.0/Reference_SD.txt", header = T)
paralogous = read.table("~/Gene_duplication/Seed_genes/paralogous_total.txt", header = T)
	


load("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Compendium_wise/development/Combined_Score.RData")
# LIST = unique(TT[TT$orthologous == 4567 & TT$orthologous.1 == 15903, "Gene.1"])
TT.copy = TT
TT = TT[TT$INTER == 4, ]
LIST = unique(TT[TT[, 14] == 1679 & TT[, 6] == 1571, "Gene"])

foreach (M = 1:length(LIST)) %dopar% {
	k = LIST[M]
	j = data.frame("name" = gene.names[neighborhood(g.grape, 1, k, mode = c("out"))[[1]],])
	j = as.character(j[,1])

	ORDER = 2
	CondRank = foreach(i = 1:length(j), .combine = rbind) %do%{
	rank(total_cond_scores[, j[i]])
	}

	Med.cond.score = apply(CondRank, 2, median)

	CondRank = rbind(c(1:dim(total_cond_scores)[1]), Med.cond.score, CondRank)
	colnames(CondRank) = colnames(data.grape)
	CondRank1 = CondRank[, order(CondRank["Med.cond.score", ], decreasing = T)]
	CondRank2 = CondRank[, order(CondRank["Med.cond.score", ], decreasing = F)]

	# CondRank1 = CondRank[, order(CondRank[ORDER + 1, ], decreasing = T)]
	# CondRank2 = CondRank[, order(CondRank[ORDER + 1, ], decreasing = F)]

	jpeg(paste("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Validation/Compendium_wise/development/Heatmap/Grape/",k,"_All_1571_1679.jpeg", sep = ""), width = 11, height = 7, units = "in", res = 500, pointsize=8)
	myImagePlot_sym(data.grape.scale[j, c(CondRank1[1, 1:10], CondRank2[1, 1:10])])
	dev.off()

	memb.comun.grape = data.frame("name" = gene.names[neighborhood(g.grape, 1, k, mode = c("out"))[[1]],])
	memb.comun.grape.para = merge(memb.comun.grape, map, by = "name")
	## Remove duplicates 
	memb.comun.grape.para = memb.comun.grape.para[!(duplicated(memb.comun.grape.para$name) 	| duplicated(memb.comun.grape.para$name, fromLast = TRUE)), ]
	memb.comun.grape.para = merge(memb.comun.grape.para, orthologous, by = "PLAZA.ID")

	j = 1571
	memb.comun.Arab = data.frame("gene_id.y" = Ref[neighborhood(g, 1, j, mode = c("out"))[[1]],1])
	# Orthologous genes between each module
	Orth1 = merge(memb.comun.Arab, memb.comun.grape.para, by = "gene_id.y")$name
	Orth1.arab = merge(memb.comun.Arab, memb.comun.grape.para, by = "gene_id.y")$gene_id.y
	
	
	# Second pair arabidopsis module to vitis module
	j = 1679
	memb.comun.Arab = data.frame("gene_id.y" = Ref[neighborhood(g, 1, j, mode = c("out"))[[1]],1])
	# Orthologous genes between each module
	Orth2 = merge(memb.comun.Arab, memb.comun.grape.para, by = "gene_id.y")$name
	Orth2.arab = merge(memb.comun.Arab, memb.comun.grape.para, by = "gene_id.y")$gene_id.y
	
	Orth = unique(intersect(as.matrix(Orth1), as.matrix(Orth2)))

	CondRank = foreach(i = 1:length(Orth1), .combine = rbind) %do%{
	rank(total_cond_scores[, Orth1[i]])
	}

	Med.cond.score = apply(CondRank, 2, median)

	CondRank = rbind(c(1:dim(total_cond_scores)[1]), Med.cond.score, CondRank)
	colnames(CondRank) = colnames(data.grape)
	CondRank1 = CondRank[, order(CondRank[ORDER + 1, ], decreasing = T)]
	CondRank2 = CondRank[, order(CondRank[ORDER + 1, ], decreasing = F)]

	jpeg(paste("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Validation/Compendium_wise/development/Heatmap/Grape/",k,"_1571.jpeg", sep = ""), width = 11, height = 7, units = "in", res = 500, pointsize=8)
	myImagePlot_sym(data.grape.scale[Orth1, c(CondRank1[1, 1:10], CondRank2[1, 1:10])])
	dev.off()


	CondRank = foreach(i = 1:length(Orth2), .combine = rbind) %do%{
	rank(total_cond_scores[, Orth2[i]])
	}

	Med.cond.score = apply(CondRank, 2, median)

	CondRank = rbind(c(1:dim(total_cond_scores)[1]), Med.cond.score, CondRank)
	colnames(CondRank) = colnames(data.grape)
	CondRank11 = CondRank[, order(CondRank[ORDER + 1, ], decreasing = T)]
	CondRank21 = CondRank[, order(CondRank[ORDER + 1, ], decreasing = F)]

	jpeg(paste("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Validation/Compendium_wise/development/Heatmap/Grape/",k,"_1679.jpeg", sep = ""), width = 11, height = 7, units = "in", res = 500, pointsize=8)
	myImagePlot_sym(data.grape.scale[Orth2, c(CondRank11[1, 1:10], CondRank21[1, 1:10])])
	dev.off()


	CondRank = foreach(i = 1:length(Orth), .combine = rbind) %do%{
	rank(total_cond_scores[, Orth[i]])
	}

	Med.cond.score = apply(CondRank, 2, median)

	CondRank = rbind(c(1:dim(total_cond_scores)[1]), Med.cond.score, CondRank)
	colnames(CondRank) = colnames(data.grape)
	CondRank111 = CondRank[, order(CondRank["Med.cond.score", ], decreasing = T)]
	CondRank211 = CondRank[, order(CondRank["Med.cond.score", ], decreasing = F)]

	jpeg(paste("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Validation/Compendium_wise/development/Heatmap/Grape/",k,"_1571_1679.jpeg", sep = ""), width = 11, height = 7, units = "in", res = 500, pointsize=8)
	myImagePlot_sym(data.grape.scale[Orth, c(CondRank111[1, 1:10], CondRank211[1, 1:10])])
	dev.off()


}


# ------------------------------------------------------------------------
# Image_scoring
# ------------------------------------------------------------------------


Image_scoring = foreach (N = 1:nrow(ORTHO), .combine = rbind) %dopar% {
	print(N)
	J = j = ORTHO$orthologous[N]
	Jpar = jpar = ORTHO$orthologous.1[N]

# foreach (N = 1:50) %dopar% {
	
	# j = Image_scoring[N, 1]
	# jpar = Image_scoring[N, 2]
	memb.comun = neighborhood(g, 1, j, mode = c("out"))[[1]]
	memb.comun.para = neighborhood(g, 1, jpar , mode = c("out"))[[1]]
	
	Inter1 = intersect(memb.comun.para, paralogous[paralogous$Ref %in% memb.comun, "Anch"])
	Inter2 = intersect(memb.comun, paralogous[paralogous$Ref %in% memb.comun.para, "Anch"])
	aa = paralogous[paralogous$Ref %in% Inter2, ]
	bb = aa[aa$Anch %in% Inter1,]
	colnames(bb)[1] = "References"
	Inter_name = merge(ref, bb, by = "References")
	colnames(Inter_name)[c(1,3)] = c("Ref", "References")
	Inter_name = merge(ref, Inter_name, by = "References")
	
	i = paste(j, jpar, sep = "_")
	write.table(Inter_name, file = paste("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Local/Compendium_wise/",NAME.spe,"/Score99/Module_Anchorpoints_V2/Module",i,".txt", sep = ""), row.names = F, quote = F, sep = "\t")


	colnames(Inter_name)[c(1,3)] = c("Ref", "References")
	Result = merge(ref, Inter_name, by = "References")
	k = paste(j, jpar, sep = "_")
	j = Result$Ref #Gene.y
	jpar = Result$References #Gene.x
	
	
	CondRank1 = foreach(i = 1:length(j), .combine = rbind) %do%{
		rank(total_cond_scores[, j[i]])
	}
	CondRank2 = foreach(i = 1:length(jpar), .combine = rbind) %do%{
		rank(total_cond_scores[, jpar[i]])
	}
	
	Med.cond.score1 = apply(CondRank1, 2, median)
	Med.cond.score2 = apply(CondRank2, 2, median)

	CondRank = rbind(c(1:dim(total_cond_scores)[1]), Med.cond.score1, CondRank1, CondRank2)
	colnames(CondRank) = colnames(data.exp.scale.med)
	CondRank0 = CondRank[,order(CondRank["Med.cond.score1", ], decreasing = T)]

	CondRank = rbind(c(1:dim(total_cond_scores)[1]), Med.cond.score2, CondRank2, CondRank1)
	colnames(CondRank) = colnames(data.exp.scale.med)
	CondRank01 = CondRank[,order(CondRank["Med.cond.score2", ], decreasing = T)]

	# jpeg(paste("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Validation/Compendium_wise/development/Heatmap/Top100/",k,".jpeg", sep = ""), width = 11, height = 7, units = "in", res = 200, pointsize=8)
	myImagePlot_sym(data.exp.scale[c(j, jpar),c(CondRank0[1,1:10],CondRank01[1,1:10])])
	# dev.off()

	A = Sig.med.zero[c(j),c(CondRank0[1,1:10])]
	A.par = Sig.med.zero[c(jpar),c(CondRank0[1,1:10])]
	Ascore = sum(A)/(nrow(A) * ncol(A))
	Ascore.par = sum(A.par)/(nrow(A) * ncol(A))

	B = Sig.med.zero[c(j), c(CondRank01[1,1:10])]
	B.par = Sig.med.zero[c(jpar), c(CondRank01[1,1:10])]
	Bscore = sum(B)/(nrow(B) * ncol(B))
	Bscore.par = sum(B.par)/(nrow(B) * ncol(B))

	c(J, Jpar, (Ascore - Ascore.par), (Bscore.par - Bscore))

}

save(Image_scoring, file = "/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Validation/Compendium_wise/abiotic/Image_scoring.RData")

load("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Validation/Compendium_wise/development/Image_scoring.RData")
Image_scoring = Image_scoring[order(Image_scoring[,3], Image_scoring[,4], decreasing = T), ]
load("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Compendium_wise/development/Combined_Score.RData")
TT["Cor.Diff"] = abs(TT$Cor1 - TT$Cor2)
TTscore = apply(TT[, -c(6:7,13:19, 25:26)], 1, sum)
TT["Score"] = TTscore
# TT = TT[order(-TT$Score),]
TT["Cor.sum"] = TT$Cor1 + TT$Cor2
TT["COR"] = TT$Cor.Diff + TT$Cor.sum



#################################################################################################################################################################################
# Gene score correlation for two example across all 14 compendia
#################################################################################################################################################################################
condition.names = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/condition_name.txt", header = T)
NAME = names(table(condition.names))

Anchorpoint = read.delim("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/Anchorpoint/Geneset1_2_53_90.txt.txt")
Anchorpoint1 = read.delim("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/Anchorpoint/AT1G50090_AT3G19710_52_90.txt")


foreach(j = 1:14) %dopar% {
	print(j)
	load(paste("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Gene_score/Compendium_wise/",as.name(NAME[j]),"/GeneScoreCorrelationMatrix.RData", sep = ""))

	AA = rbind(cbind(COR[as.character(Anchorpoint$Pollen), as.character(Anchorpoint$Pollen)], COR[as.character(Anchorpoint$Pollen), as.character(Anchorpoint$Leaf)]),
	+ cbind(COR[as.character(Anchorpoint$Leaf), as.character(Anchorpoint$Pollen)], COR[as.character(Anchorpoint$Leaf), as.character(Anchorpoint$Leaf)]))
	BB = rbind(cbind(COR[as.character(Anchorpoint1$Root), as.character(Anchorpoint1$Root)], COR[as.character(Anchorpoint1$Root), as.character(Anchorpoint1$Leaf)]),
	+ cbind(COR[as.character(Anchorpoint1$Leaf), as.character(Anchorpoint1$Root)], COR[as.character(Anchorpoint1$Leaf), as.character(Anchorpoint1$Leaf)]))
	
	write.table(AA, file = paste("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Gene_score/Compendium_wise/",as.name(NAME[j]),"/GeneScoreCorrelation_Pollen_Leaf.txt", sep = ""),
		quote = F, sep = "\t")
	write.table(BB, file = paste("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Gene_score/Compendium_wise/",as.name(NAME[j]),"/GeneScoreCorrelation_Root_Leaf.txt", sep = ""),
		quote = F, sep = "\t")
	
	jpeg(paste("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Gene_score/Compendium_wise/",as.name(NAME[j]),"/GeneScoreCorrelation_Pollen_Leaf.jpeg", sep = ""), width = 11, height = 7, units = "in", res = 500, pointsize=8)
	myImagePlot(AA)
	dev.off()

	jpeg(paste("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Gene_score/Compendium_wise/",as.name(NAME[j]),"/GeneScoreCorrelation_Root_Leaf.jpeg", sep = ""), width = 11, height = 7, units = "in", res = 500, pointsize=8)
	myImagePlot(BB)
	dev.off()
}
#################################################################################################################################################################################
