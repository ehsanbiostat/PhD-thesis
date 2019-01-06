
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

data.exp = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/Integrated_SD.txt", colClasses = "numeric", nrows = 19285, comment.char = "", header = T)
data.exp.scale = data.exp[, condition.names == NAME.spe]
data.exp.scale = apply(data.exp.scale,1,scale)
data.exp.scale = t(data.exp.scale)
data.exp.scale = apply(data.exp.scale, 2, scale)



data.exp.scale.med = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/Integrated_SD_Norm_Med.txt", colClasses = "numeric", nrows = 19285, comment.char = "", header = T)

data.exp.scale.med = data.exp.scale.med[, condition.names == NAME.spe]

Sig = ifelse(data.exp.scale > 0, 1, ifelse(data.exp.scale < 0, -1, 0))
Sig.med = ifelse(data.exp.scale.med > 0, 1, ifelse(data.exp.scale.med < 0, -1, 0))

Sig.med.zero = ifelse(data.exp.scale.med > 0, 1, ifelse(data.exp.scale.med < 0, 0, 0))


ref = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/Reference_SD.txt", header = T)
# paralogous = read.table("~/Gene_duplication/Seed_genes/paralogous_total.txt", header = T)
paralogous = read.table("~/Gene_duplication/Seed_genes/paralogous_NoFilter.txt", header = T)
## Load generated graph from SA result
load(paste("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Gene_score/Compendium_wise/",NAME.spe,"/Graph_V2.RData", sep = ""))


ORTHO = read.table("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Validation/Compendium_wise/development/4.Para.list.txt", header = T)


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
data.grape = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/Vitis_vinifera/GE.mean.SD.txt", nrows = 19589, comment.char = "", header = T)
data.grape.scale = apply(data.grape, 1, scale)
data.grape.scale = t(data.grape.scale)
colnames(data.grape.scale) = colnames(data.grape)

load("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Graph_V2.R")
	g.grape = g

gene.names = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/Vitis_vinifera/gene.names.txt")
cond.devel = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/Vitis_vinifera/development.annotaion.txt")

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



list = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Scripts/General/order.gene_score.Grape.txt")
list = as.vector(unlist(list))
## Read text or RData file with read.table or load function
setwd("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Raw/Gene_score/")
# C = lapply(list, read.table, sep= "\t" , header=T)
G = lapply(list, read.table, header = T)
total_gene_scores = matrix(NA, length(unlist(G[1])), length(list))
for (i in 1:length(list)){
	print(i)
	total_gene_scores[, i] = as.matrix(unlist(G[i]))
}
colnames(total_gene_scores) = rownames(total_gene_scores) = gene.names[, 1]


orthologous = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Orthology/Grape-Arabidopsis/All.txt", header = T)
gene.names = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/Vitis_vinifera/gene.names.txt")
map = read.delim("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Orthology/Grape-Arabidopsis/Map.12x.plaza.txt")
Ref = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/Reference_SD.txt", header = T)
paralogous = read.table("~/Gene_duplication/Seed_genes/paralogous_total.txt", header = T)
	



LIST = unique(TT[TT$orthologous == 1679 & TT$orthologous.1 == 10809, "Gene.1"])

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

	CondRank1 = CondRank[, order(CondRank[ORDER + 1, ], decreasing = T)]
	CondRank2 = CondRank[, order(CondRank[ORDER + 1, ], decreasing = F)]
	sum(data.grape.scale[rownames(data.grape.scale) %in% as.character(RootLeaf$name), c(CondRank1[1, 1:10])] > 0)/sum(data.grape.scale[rownames(data.grape.scale) %in% as.character(RootLeaf$name), c(CondRank2[1, 1:10])] > 0)

	jpeg(paste("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Validation/Compendium_wise/development/Heatmap/Grape/",k,"_1679_10809.jpeg", sep = ""), width = 11, height = 7, units = "in", res = 500, pointsize=8)
	myImagePlot_sym(data.grape.scale[j, c(CondRank1[1, 1:10], CondRank2[1, 1:10])])
	dev.off()

	memb.comun.grape = data.frame("name" = gene.names[neighborhood(g.grape, 1, k, mode = c("out"))[[1]],])
	memb.comun.grape.para = merge(memb.comun.grape, map, by = "name")
	## Remove duplicates 
	memb.comun.grape.para = memb.comun.grape.para[!(duplicated(memb.comun.grape.para$name) 	| duplicated(memb.comun.grape.para$name, fromLast = TRUE)), ]
	memb.comun.grape.para = merge(memb.comun.grape.para, orthologous, by = "PLAZA.ID")

	j = 1679
	memb.comun.Arab = data.frame("gene_id.y" = Ref[neighborhood(g, 1, j, mode = c("out"))[[1]],1])
	# Orthologous genes between each module
	Orth1 = merge(memb.comun.Arab, memb.comun.grape.para, by = "gene_id.y")$name
	Orth1.arab = merge(memb.comun.Arab, memb.comun.grape.para, by = "gene_id.y")$gene_id.y
	
	
	# Second pair arabidopsis module to vitis module
	j = 10809
	memb.comun.Arab = data.frame("gene_id.y" = Ref[neighborhood(g, 1, j, mode = c("out"))[[1]],1])
	# Orthologous genes between each module
	Orth2 = merge(memb.comun.Arab, memb.comun.grape.para, by = "gene_id.y")$name
	Orth2.arab = merge(memb.comun.Arab, memb.comun.grape.para, by = "gene_id.y")$gene_id.y
	
	Orth = unique(c(as.matrix(Orth1), as.matrix(Orth2)))

	CondRank = foreach(i = 1:length(Orth1), .combine = rbind) %do%{
	rank(total_cond_scores[, Orth1[i]])
	}

	Med.cond.score = apply(CondRank, 2, median)

	CondRank = rbind(c(1:dim(total_cond_scores)[1]), Med.cond.score, CondRank)
	colnames(CondRank) = colnames(data.grape)
	CondRank1 = CondRank[, order(CondRank[ORDER + 1, ], decreasing = T)]
	CondRank2 = CondRank[, order(CondRank[ORDER + 1, ], decreasing = F)]

	jpeg(paste("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Validation/Compendium_wise/development/Heatmap/Grape/",k,"_1679.jpeg", sep = ""), width = 11, height = 7, units = "in", res = 500, pointsize=8)
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

	jpeg(paste("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Validation/Compendium_wise/development/Heatmap/Grape/",k,"_10809.jpeg", sep = ""), width = 11, height = 7, units = "in", res = 500, pointsize=8)
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

	jpeg(paste("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Validation/Compendium_wise/development/Heatmap/Grape/",k,"_1679_10809.jpeg", sep = ""), width = 11, height = 7, units = "in", res = 500, pointsize=8)
	myImagePlot_sym(data.grape.scale[Orth, c(CondRank111[1, 1:10], CondRank211[1, 1:10])])
	dev.off()

}
	


# ------------------------------------------------------------------------
# Image_scoring
# ------------------------------------------------------------------------

ORTHO = read.table("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Validation/Compendium_wise/development/5.Para.list.txt", header = T)
Image_scoring = foreach (N = 1:nrow(ORTHO), .combine = rbind) %dopar% {
	print(N)
	J = j = ORTHO$orthologous[N]
	Jpar = jpar = ORTHO$orthologous.1[N]

	# foreach (N = 1:nrow(Image_scoring_All)) %dopar% {
	
	#  	j = Image_scoring_All[N, 1]
	#  	jpar = Image_scoring_All[N, 2]
	  	memb.comun = neighborhood(g,1,j, mode = c("out"))[[1]]
		memb.comun.para = neighborhood(g,1,jpar, mode = c("out"))[[1]]
		
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

	#  	jpeg(paste("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Validation/Compendium_wise/development/Heatmap/TopAll/",k,".jpeg", sep = ""), width = 11, height = 7, units = "in", res = 200, pointsize=8)
	#  	myImagePlot_sym(data.exp.scale[c(j, jpar),c(CondRank0[1,1:10],CondRank01[1,1:10])])
	#  	dev.off()
	# }

	A = Sig.med.zero[c(j),c(CondRank0[1,1:10])]
	A.par = Sig.med.zero[c(jpar),c(CondRank0[1,1:10])]
	Ascore = sum(A)/(nrow(A) * ncol(A))
	Ascore.par = sum(A.par)/(nrow(A) * ncol(A))

	B = Sig.med.zero[c(j), c(CondRank01[1,1:10])]
	B.par = Sig.med.zero[c(jpar), c(CondRank01[1,1:10])]
	Bscore = sum(B)/(nrow(B) * ncol(B))
	Bscore.par = sum(B.par)/(nrow(B) * ncol(B))

	c(J, Jpar, (Ascore - Ascore.par), (Bscore.par - Bscore), ORTHO$Paralogous.cluster[N])

}

save(Image_scoring, file = "/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Validation/Compendium_wise/development/Image_scoring.RData")

load("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Validation/Compendium_wise/development/Image_scoring.RData")
Image_scoring = Image_scoring[order(Image_scoring[,5], Image_scoring[,3], Image_scoring[,4], decreasing = T), ]

Image_scoring_All = rbind(matrix(Image_scoring[Image_scoring[,5] == 9,][1,], 1), 
	Image_scoring[Image_scoring[,5] == 8,][1:3,], 
	Image_scoring[Image_scoring[,5] == 7,][1:6,], 
	Image_scoring[Image_scoring[,5] == 6,][1:12,], 
	Image_scoring[Image_scoring[,5] == 5,][1:24,], 
	Image_scoring[Image_scoring[,5] == 4,][1:48,]
)



# ------------------------------------------------------------------------
# Image_scoring for anchorpoints
# ------------------------------------------------------------------------

ORTHO = read.table("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Validation/Compendium_wise/development/Anch.list.txt", header = T)

ORTHO.anch = foreach(i = 1:nrow(paralogous), .combine = rbind) %dopar% {
	print(i)
	ORTHO[ORTHO[, 1]%in% paralogous$Ref[i] & ORTHO[, 2] %in% paralogous$Anch[i] ,]
}

top = read.table("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Validation/Compendium_wise/development/Heatmap/TOP.txt")
top = read.delim(paste("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/",NAME[k],"/3threshold_20Jaccard.txt", sep = ""))
ORTHO = top[top$Threshold == 95,]

# Threshold 95%
Score = read.delim("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Local/Compendium_wise/development/Score95/Anch/Scores_V2/Anch_Score.txt")
Score = Score[Score$Jaccard == 0,]
Score = Score[order(-Score$Normalized.Paralogous.cluster), ]

Image_scoring = foreach (N = 1:nrow(ORTHO), .combine = rbind) %dopar% {
	print(N)
	J = j = ORTHO$Gene[N]
	Jpar = jpar = ORTHO$Paralogous[N]
  	memb.comun = neighborhood(g,1,j, mode = c("out"))[[1]]
	memb.comun.para = neighborhood(g,1,jpar, mode = c("out"))[[1]]
		
	Inter1 = intersect(memb.comun.para, paralogous[paralogous$Ref %in% memb.comun, "Anch"])
	Inter2 = intersect(memb.comun, paralogous[paralogous$Ref %in% memb.comun.para, "Anch"])
	aa = paralogous[paralogous$Ref %in% Inter2, ]
	bb = aa[aa$Anch %in% Inter1,]
	colnames(bb)[1] = "References"
	Inter_name = merge(ref, bb, by = "References")
	colnames(Inter_name)[c(1,3)] = c("Ref", "References")
	Inter_name = merge(ref, Inter_name, by = "References")
		
	colnames(Inter_name)[c(1,3)] = c("Ref", "References")
	Result = merge(ref, Inter_name, by = "References")
	# k = paste(ref[j, 1], ref[jpar, 1], Image_scoring[N, 5], sep = "_")
	j = Result$References #Gene.y
	jpar = Result$Ref #Gene.x
		
		
	CondRank1 = foreach(i = 1:length(memb.comun), .combine = rbind) %do%{
		rank(total_cond_scores[, memb.comun[i]])
	}
	CondRank2 = foreach(i = 1:length(memb.comun.para), .combine = rbind) %do%{
		rank(total_cond_scores[, memb.comun.para[i]])
	}
	
	# WW = foreach(JJ = which(memb.comun %in% Geneset1), .combine =rbind) %do% {
		# W = foreach(KK = which(memb.comun.para %in% Geneset2), .combine =rbind) %do% {
	# ORDER = JJ 
	ORDER = 1
	# Med.cond.score1 = apply(CondRank1, 2, median)
	# Med.cond.score2 = apply(CondRank2, 2, median)

		
	CondRank = rbind(c(1:dim(total_cond_scores)[1]), CondRank1, CondRank2)
	colnames(CondRank) = colnames(data.exp.scale.med)
	CondRank0 = CondRank[,order(CondRank[ORDER + 1, ], decreasing = T)]

	# ORDER = KK
	CondRank = rbind(c(1:dim(total_cond_scores)[1]), CondRank2, CondRank1)
	colnames(CondRank) = colnames(data.exp.scale.med)
	CondRank01 = CondRank[,order(CondRank[ORDER + 1, ], decreasing = T)]
		
	A = Sig[c(j),c(CondRank0[1,1:10])]
	A.par = Sig[c(jpar),c(CondRank0[1,1:10])]
	Ascore = sum(A - A.par)/(nrow(A) * ncol(A) * 2)
	if(length(Ascore) == 0) Ascore = NA
	B = Sig[c(j), c(CondRank01[1,1:10])]
	B.par = Sig[c(jpar), c(CondRank01[1,1:10])]
	Bscore = sum(B.par - B)/(nrow(B) * ncol(B) * 2)
	if(length(Bscore) == 0) Bscore = NA
	c(J, Jpar, Ascore, Bscore, ORTHO$Paralogous.cluster[N])
	# c(JJ, KK, Ascore, Bscore)
	# }
	# W
	# }

}

Image_scoring = cbind(Image_scoring, (Image_scoring[,3] + Image_scoring[,4]))
Image_scoring = Image_scoring[order(Image_scoring[,6], Image_scoring[, 3], Image_scoring[,4], decreasing = T), ]

leaf = read.table("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/Module/Leaf_90.txt")
flower = read.table("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/Module/Pollen_90.txt")


leaf = read.table("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/Module/AT3G19710_Leaf_90.txt")
Root = read.table("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/Module/AT1G50090_Root_90.txt")
rownames(ref) = ref[,1]

foreach (N = 1:nrow(Image_scoring)) %dopar% {
	
	j = Image_scoring[N, 1]
	jpar = Image_scoring[N, 2]


	# Annot = foreach (N = 1:nrow(top), .combine = rbind) %dopar% {
	# j = top[N, 1]
	# jpar = top[N, 2]
	  	memb.comun = neighborhood(g,1,j, mode = c("out"))[[1]]
		memb.comun.para = neighborhood(g,1,jpar, mode = c("out"))[[1]]

		memb.comun = ref[as.character(leaf[,1]),2]
		memb.comun.para = ref[as.character(Root[,1]),2]
		
		Inter1 = intersect(memb.comun.para, paralogous[paralogous$Ref %in% memb.comun, "Anch"])
		Inter2 = intersect(memb.comun, paralogous[paralogous$Ref %in% memb.comun.para, "Anch"])
		aa = paralogous[paralogous$Ref %in% Inter2, ]
		bb = aa[aa$Anch %in% Inter1,]
		colnames(bb)[1] = "References"
		Inter_name = merge(ref, bb, by = "References")
		colnames(Inter_name)[c(1,3)] = c("Ref", "References")
		Inter_name = merge(ref, Inter_name, by = "References")
		
		i = paste(ref[j, 1], ref[jpar, 1], Image_scoring[N, 5], Image_scoring[N, 6], Image_scoring$Threshold[N], sep = "_")
		# i = paste(ref[j, 1], ref[jpar, 1], top[N, 5], top[N, 6], sep = "_")
		# write.table(Inter_name, file = paste("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Local/Compendium_wise/",NAME.spe,"/Score99/Module_Anchorpoints_V2/Module",i,".txt", sep = ""), row.names = F, quote = F, sep = "\t")
		# write.table(Inter_name, file = paste("/group/biocomp/users/ehsab/Gene_duplication/Results/SA/Local/Validation/Compendium_wise/Selected_V4/development/95/Genes/Anchorpoint/",i,".txt", sep = ""), row.names = F, quote = F, sep = "\t")
		# write.table(ref[memb.comun, 1], file = paste("/group/biocomp/users/ehsab/Gene_duplication/Results/SA/Local/Validation/Compendium_wise/Selected_V4/development/95/Genes/Module/",ref[j, 1],".txt", sep = ""), col.names = F, row.names = F, quote = F, sep = "\t")
		# write.table(ref[memb.comun.para, 1], file = paste("/group/biocomp/users/ehsab/Gene_duplication/Results/SA/Local/Validation/Compendium_wise/Selected_V4/development/95/Genes/Module/",ref[jpar, 1],".txt", sep = ""), col.names = F, row.names = F, quote = F, sep = "\t")

		# write.table(Inter_name, file = paste("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/Anchorpoint/",i,".txt", sep = ""), row.names = F, quote = F, sep = "\t")
		# write.table(ref[memb.comun, 1], file = paste("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/Module/",ref[j, 1],"_",Image_scoring$Threshold[N],".txt", sep = ""), col.names = F, row.names = F, quote = F, sep = "\t")
		# write.table(ref[memb.comun.para, 1], file = paste("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/Module/",ref[jpar, 1],"_",Image_scoring$Threshold[N],".txt", sep = ""), col.names = F, row.names = F, quote = F, sep = "\t")


		colnames(Inter_name)[c(1,3)] = c("Ref", "References")
		Result = merge(ref, Inter_name, by = "References")
		# k = paste(ref[j, 1], ref[jpar, 1], Image_scoring[N, 5], sep = "_")
		# k = paste(ref[j, 1], ref[jpar, 1], top[N, 5], top[N, 6], sep = "_")
		k = paste(ref[j, 1], ref[jpar, 1], Image_scoring[N, 5], Image_scoring[N, 6], Image_scoring$Threshold[N], sep = "_")
		j = Result$References #Gene.y
		jpar = Result$Ref #Gene.x
		
		
		CondRank1 = foreach(i = 1:length(memb.comun), .combine = rbind) %do%{
			rank(total_cond_scores[, memb.comun[i]])
		}
		CondRank2 = foreach(i = 1:length(memb.comun.para), .combine = rbind) %do%{
			rank(total_cond_scores[, memb.comun.para[i]])
		}
		
		ORDER = 1 
		# Med.cond.score1 = apply(CondRank1[which(memb.comun %in% Geneset1),], 2, median)
		# Med.cond.score2 = apply(CondRank2[which(memb.comun.para %in% Geneset2),], 2, median)

		# CondRank = rbind(c(1:length(CondRank1)), data.frame(t(CondRank1)), data.frame(t(CondRank2)))
		CondRank = rbind(c(1:dim(total_cond_scores)[1]), CondRank1, CondRank2)
		colnames(CondRank) = colnames(data.exp.scale.med)
		CondRank0 = CondRank[,order(CondRank[ORDER + 1, ], decreasing = T)]
		# CondRank0 = CondRank[,order(Med.cond.score1, decreasing = T)]

		CondRank = rbind(c(1:dim(total_cond_scores)[1]), CondRank2, CondRank1)
		colnames(CondRank) = colnames(data.exp.scale.med)
		CondRank01 = CondRank[,order(CondRank[ORDER + 1, ], decreasing = T)]
		# CondRank01 = CondRank[,order(Med.cond.score2, decreasing = T)]
		
	 	
		# Entire module heatmap
		# jpeg(paste("/group/biocomp/users/ehsab/Gene_duplication/Results/SA/Local/Validation/Compendium_wise/Selected_V4/development/95/Plot/Module/Top10/",k,".jpeg", sep = ""), width = 11, height = 7, units = "in", res = 200, pointsize=8)
		jpeg(paste("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Plot/Module/Top10/",k,".jpeg", sep = ""), width = 11, height = 7, units = "in", res = 200, pointsize=8)
	   	myImagePlot_sym(data.exp.scale[c(memb.comun, memb.comun.para),c(CondRank0[1,1:10],CondRank01[1,1:10])])
	   	dev.off()
		


	 	# Anchorpoints module heatmap
	 	# jpeg(paste("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Validation/Compendium_wise/development/Heatmap/Anchorpoints/",k,".jpeg", sep = ""), width = 11, height = 7, units = "in", res = 200, pointsize=8)
	 	# jpeg(paste("/group/biocomp/users/ehsab/Gene_duplication/Results/SA/Local/Validation/Compendium_wise/Selected_V4/development/95/Plot/Anchorpoint/Top10/",k,".jpeg", sep = ""), width = 11, height = 7, units = "in", res = 200, pointsize=8)
	 	jpeg(paste("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Plot/Anchorpoint/Top10/",k,".jpeg", sep = ""), width = 11, height = 7, units = "in", res = 200, pointsize=8)
	 	myImagePlot_sym(data.exp.scale[c(j, jpar),c(CondRank0[1,1:10],CondRank01[1,1:10])])
	 	dev.off()
		
	 	# ------------------------------------------------------------------
	 	# CORNET annotation
	 	k = 134
		Annot.j = foreach(i = 1:k, .combine = c) %do% {
			strsplit(substring(colnames(CondRank0)[1:k], 2), "[.]")[[i]][1]
		}

		Annot.jpar = foreach(i = 1:k, .combine = c) %do% {
			strsplit(substring(colnames(CondRank01)[1:k], 2), "[.]")[[i]][1]
		}
		# ------------------------------------------------------------------

		# data.frame(cbind(Annot.j, Annot.jpar), memb.comun[1], memb.comun.para[1])
	}
}

load("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/PO_CORNET_Annotation.RData")
COR.annot.PO[names(COR.annot.PO) %in% Annot.j]


save(Image_scoring, file = "/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Validation/Compendium_wise/development/Image_scoring.RData")
save(Annot, file = "/scratch/iwt/ehsab/Gene_duplication/Results/SA/Local/Compendium_wise/development/CORNET.annot.RData")










# CORNET annotation
CORNET.annot = read.delim("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/cornet_allMA_desc_120509.txt")
which(substring(strsplit(as.character(CORNET.annot$Description[2]), ",")[[1]], 1, 3) == " PO")
CORNET.annot[CORNET.annot$Name %in% Annot.j,]



AnchorpointShared = foreach (N = 1:nrow(Image_scoring), .combine = rbind) %dopar% {
	
	j = Image_scoring[N, 1]
	jpar = Image_scoring[N, 2]


	  	memb.comun = neighborhood(g,1,j, mode = c("out"))[[1]]
		memb.comun.para = neighborhood(g,1,jpar, mode = c("out"))[[1]]
		
		Inter1 = intersect(memb.comun.para, paralogous[paralogous$Ref %in% memb.comun, "Anch"])
		Inter2 = intersect(memb.comun, paralogous[paralogous$Ref %in% memb.comun.para, "Anch"])
		aa = paralogous[paralogous$Ref %in% Inter2, ]
		bb = aa[aa$Anch %in% Inter1,]
		colnames(bb)[1] = "References"
		Inter_name = merge(ref, bb, by = "References")
		colnames(Inter_name)[c(1,3)] = c("Ref", "References")
		Inter_name = merge(ref, Inter_name, by = "References")
		i = paste(j, jpar, Image_scoring[N, 5], 95, sep = "_")
		Inter_name["Group"] = i
		Inter_name
}
		


AnnotCornet = foreach (N = 1:nrow(Image_scoring), .combine = rbind) %dopar% {
	
	j = Image_scoring[N, 1]
	jpar = Image_scoring[N, 2]


	  	memb.comun = neighborhood(g,1,j, mode = c("out"))[[1]]
		memb.comun.para = neighborhood(g,1,jpar, mode = c("out"))[[1]]
		
		Inter1 = intersect(memb.comun.para, paralogous[paralogous$Ref %in% memb.comun, "Anch"])
		Inter2 = intersect(memb.comun, paralogous[paralogous$Ref %in% memb.comun.para, "Anch"])
		aa = paralogous[paralogous$Ref %in% Inter2, ]
		bb = aa[aa$Anch %in% Inter1,]
		colnames(bb)[1] = "References"
		Inter_name = merge(ref, bb, by = "References")
		colnames(Inter_name)[c(1,3)] = c("Ref", "References")
		Inter_name = merge(ref, Inter_name, by = "References")		
		colnames(Inter_name)[c(1,3)] = c("Ref", "References")
		Result = merge(ref, Inter_name, by = "References")
		k = paste(ref[j, 1], ref[jpar, 1], Image_scoring[N, 5], Image_scoring[N, 6], Image_scoring$Threshold[N], sep = "_")
		j = Result$References #Gene.y
		jpar = Result$Ref #Gene.x
		
		
		CondRank1 = foreach(i = 1:length(memb.comun), .combine = rbind) %do%{
			rank(total_cond_scores[, memb.comun[i]])
		}
		CondRank2 = foreach(i = 1:length(memb.comun.para), .combine = rbind) %do%{
			rank(total_cond_scores[, memb.comun.para[i]])
		}
		
		ORDER = 1 
		# Med.cond.score1 = apply(CondRank1, 2, median)
		# Med.cond.score2 = apply(CondRank2, 2, median)

		
		CondRank = rbind(c(1:dim(total_cond_scores)[1]), CondRank1, CondRank2)
		colnames(CondRank) = colnames(data.exp.scale.med)
		CondRank0 = CondRank[,order(CondRank[ORDER + 1, ], decreasing = T)]

		CondRank = rbind(c(1:dim(total_cond_scores)[1]), CondRank2, CondRank1)
		colnames(CondRank) = colnames(data.exp.scale.med)
		CondRank01 = CondRank[,order(CondRank[ORDER + 1, ], decreasing = T)]
		
	 	# ------------------------------------------------------------------
	 	# CORNET annotation
	 	k = 10
		Annot.j = foreach(i = 1:k, .combine = c) %do% {
			strsplit(substring(colnames(CondRank0)[1:k], 2), "[.]")[[i]][1]
		}

		Annot.jpar = foreach(i = 1:k, .combine = c) %do% {
			strsplit(substring(colnames(CondRank01)[1:k], 2), "[.]")[[i]][1]
		}
		# ------------------------------------------------------------------

		data.frame(rbind(Annot.j, Annot.jpar), "Seed" = rbind(memb.comun[1], memb.comun.para[1]))

}

Jaccard.Cornet = foreach (i = 1:nrow(AnnotCornet), .combine = rbind) %do% {
	JJ = foreach(j = 1:nrow(AnnotCornet), .combine = c) %do% {
		length(intersect(as.vector(unlist(AnnotCornet[i,1:10])), as.vector(unlist(AnnotCornet[j,1:10]))))/length(union(as.vector(unlist(AnnotCornet[i,1:10])), as.vector(unlist(AnnotCornet[j,1:10]))))
	}
	JJ
}













# Two examples
rownames(data.exp.scale) = ref[,1]
colnames(data.exp.scale) = colnames(data.exp.scale.med)
load("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Condition_score/Compendium_wise/development/Raw/Geneset1_leaf.RData")
cond_scoresLeaf = cond_scores
load("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Condition_score/Compendium_wise/development/Raw/Geneset1_pollen.RData")
cond_scoresPollen = cond_scores

Leaf = read.table("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/Module/Leaf_90.txt")

Anch = read.delim("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/Anchorpoint/Geneset1_2_53_90.txt.txt")

jpeg("/group/biocomp/users/ehsab/Gene_duplication/Results/SA/Local/Validation/Compendium_wise/Selected_V5/2candidates/Heatmap/Leaf_Pollen_90_Anch1.jpeg", width = 11, height = 7, units = "in", res = 500, pointsize=8)
myImagePlot_sym(data.exp.scale[c(as.character(Anch$Leaf), as.character(Anch$Pollen)) , c(names(cond_scoresLeaf[order(-cond_scoresLeaf[,1]),])[1:10], names(cond_scoresPollen[order(-cond_scoresPollen[,1]),])[1:10])])
dev.off()

AA = cbind(data.exp.scale[c(as.character(Anch$Leaf), as.character(Anch$Pollen)) , c(names(cond_scoresLeaf[order(-cond_scoresLeaf[,1]),])[1:10], names(cond_scoresPollen[order(-cond_scoresPollen[,1]),])[1:10])], c(stability.AnchorRoot, stability.AnchorRootLeaf))


LeafTop = data.exp.scale[c(as.character(Anch$Leaf)), names(cond_scoresLeaf[order(-cond_scoresLeaf[,1]),])[1:10]]

Anch = Anch[order(Anch$WGDevent),]
edgelist = as.matrix(Anch[,c(2,4)])
Vertices = unique(c(as.character(Anch$Leaf), as.character(Anch$Pollen)))
VV1 = c(Vertices[1:52], Vertices[length(Vertices):53])
VV = unique(VV1)
edgeColor = rep("green2", nrow(edgelist))
edgeColor[which(Anch$WGDevent == "recent")] = "red"
above = rep(TRUE, nrow(edgelist))
above[sample(1:length(above), 53, replace = F)] = "FALSE"


source("/ngsprojects/iwt/ehsab/Gene_duplication/Scripts/arcplot/arcplot.r")

jpeg("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/ARC.jpeg", width=8, height=3, units="in", res=800, pointsize = 5)
arcplot(edgelist, vertices = VV, col.arcs = edgeColor, col.nodes = "black", las = 3, lwd.arc = 1, lwd.nodes= 2, line = 1, lty.arcs = 1)
dev.off()




3023
FF = foreach(i =4:3023, .combine = rbind) %dopar% {
	print(i)
	load(paste("Threshold",i,".RData", sep = ""))
	AA = data.frame("Mean" = apply(Result, 1, mean), "SD" = apply(Result, 1, sd))
	Anc = data.frame("Anc" = AA[seq.int(1, nrow(AA), 2),], "Jac" = AA[seq.int(2, nrow(AA), 2),], "Module1" = i, "Module2" = 4:3026)
	write.table(Anc, file = paste("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Local/Compendium_wise/Randomization/MeanSD/",i,".txt", sep = ""), sep = "\t", row.names = F, quote = F)
	Anc
}


#####################################################################################################################################################################################
# Stability ang gene scores with heatmap
#
#####################################################################################################################################################################################

Pollen = read.delim("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/StabilityPollenLeaf.txt")
Root = read.delim("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/StabilityRootLeaf.txt")
PollenRoot = Pollen[Pollen$Ref %in% Root$Ref & Pollen$References %in% Root$References,]
RootPollen = Root[Root$Ref %in% Pollen$Ref & Root$References %in% Pollen$References,]
PollenRoot = PollenRoot[order(PollenRoot$Ref, PollenRoot$References),]
RootPollen = RootPollen[order(RootPollen$Ref, RootPollen$References),]
PollenRootTotal = data.frame("Root/Pollen" = RootPollen[, 2],"Leaf" = RootPollen[, 4],
	"stabilityRoot" = RootPollen[, 7], "stabilityPollen" = PollenRoot[,c(7)], "stabilityLeafRoot" = RootPollen[, 8],
	"stabilityLeafPollen" = PollenRoot[,8], "RootGeneScore" = RootPollen[,9], "PollenGeneScore" = PollenRoot[,9],
	"LeafRootGeneScore" = RootPollen[,10], "LeafPollenGeneScore" = PollenRoot[,10], "RootGeneScoreRank" = RootPollen[,11],
	"PollenGeneScoreRank" = PollenRoot[,11], "LeafRootGeneScoreRank" = RootPollen[,12], "LeafPollenGeneScoreRank" = PollenRoot[,12])
write.table(PollenRootTotal, file = "/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/StabilityPollenRoot.txt",
	row.names = F, quote = F, sep = "\t")


# I checked them manually and in some cases the conclusions are not the same
PollenRootTotal["Module"] = ifelse(PollenRootTotal$stabilityRoot > PollenRootTotal$stabilityPollen & PollenRootTotal$stabilityLeafRoot > PollenRootTotal$stabilityLeafPollen, "Root",
	ifelse(PollenRootTotal$stabilityRoot < PollenRootTotal$stabilityPollen & PollenRootTotal$stabilityLeafRoot < PollenRootTotal$stabilityLeafPollen, "Pollen",
		ifelse(PollenRootTotal$stabilityRoot > PollenRootTotal$stabilityPollen & PollenRootTotal$stabilityLeafRoot < PollenRootTotal$stabilityLeafPollen & PollenRootTotal$LeafRootGeneScoreRank < PollenRootTotal$LeafPollenGeneScoreRank, "Root",
			ifelse(PollenRootTotal$stabilityRoot < PollenRootTotal$stabilityPollen & PollenRootTotal$stabilityLeafRoot > PollenRootTotal$stabilityLeafPollen & PollenRootTotal$LeafRootGeneScoreRank > PollenRootTotal$LeafPollenGeneScoreRank, "Pollen", "Pollen/Root"))))



PollenStability = c(ifelse(Pollen$stabilityLeaf < 0.4, -4, 4), ifelse(Pollen$stabilityPollen < 0.4, -4, 4))
RootStability = c(ifelse(Root$stabilityLeaf < 0.3, -4, 4), ifelse(Root$stabilityRoot < 0.3, -4, 4))

load("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Gene_score/Compendium_wise/development/Raw/G2858.RData")
RootS = gene_scores
RootSR = rank(-1 * RootS)
load("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Gene_score/Compendium_wise/development/Raw/G9603.RData")
LeafS = gene_scores
LeafSR = rank(-1 * LeafS)
Root["RootGeneScoreRank"] = RootSR[Root$Ref]
Root["LeafGeneScoreRank"] = LeafSR[Root$References]
RootGeneScore = c(ifelse(LeafS[Root$References,] < 0.01, -4, 4), ifelse(RootS[Root$Ref,] < 0.01, -4, 4))


load("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Gene_score/Compendium_wise/development/Raw/Geneset2.RData")
PollenS = gene_scores
PollenSR = rank(-1 * PollenS)
load("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Gene_score/Compendium_wise/development/Raw/Geneset1.RData")
LeafS = gene_scores
LeafSR = rank(-1 * LeafS)
Pollen["PollenGeneScoreRank"] = PollenSR[Pollen$Ref]
Pollen["LeafGeneScoreRank"] = LeafSR[Pollen$References]
PollenGeneScore = c(ifelse(LeafS[Pollen$References,] < 0.01, -4, 4), ifelse(PollenS[Pollen$Ref,] < 0.01, -4, 4))



RootHeatmap = cbind(data.exp.scale[c(j, jpar),c(CondRank0[1,1:10],CondRank01[1,1:10])], RootStability, RootGeneScore)
PollenHeatmap = cbind(data.exp.scale[c(as.character(Anch$Leaf), as.character(Anch$Pollen)),
	c(names(cond_scoresLeaf[order(-cond_scoresLeaf[,1]),])[1:10], names(cond_scoresPollen[order(-cond_scoresPollen[,1]),])[1:10])], PollenStability, PollenGeneScore)
colnames(PollenHeatmap)[c(21, 22)] = colnames(RootHeatmap)[c(21, 22)] = c("Stability", "Gene score")


jpeg("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Plot/Anchorpoint/Top10/StabilityPollenLeaf.jpeg", width=11, height=7, units="in", res=800, pointsize = 7)
myImagePlot_sym(PollenHeatmap)
dev.off()

jpeg("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Plot/Anchorpoint/Top10/StabilityRootLeaf.jpeg", width=11, height=7, units="in", res=800, pointsize = 7)
myImagePlot_sym(RootHeatmap)
dev.off()
