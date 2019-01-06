

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

Anchorpoint = read.delim("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/Anchorpoint/Geneset1_2_53_90.txt.txt")
Anchorpoint1 = read.delim("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/Anchorpoint/AT1G50090_AT3G19710_52_90.txt")

colnames(Anchorpoint)[2] = "Root/Pollen"
colnames(Anchorpoint1)[4] = "Root/Pollen"
Anch = rbind(Anchorpoint, Anchorpoint1)
Anch = Anch[!(duplicated(Anch)), ]


memb.comun = Anchorpoint1$Ref
memb.comun.para = Anchorpoint1$References

memb.comun = Anchorpoint$Ref
memb.comun.para = Anchorpoint$References

library(NMF)

function(memb.comun, memb.comun.para) {
	
	CondRank1 = foreach(i = 1:length(memb.comun), .combine = rbind) %do%{
		rank(total_cond_scores[, memb.comun[i]])
	}
	CondRank2 = foreach(i = 1:length(memb.comun.para), .combine = rbind) %do%{
		rank(total_cond_scores[, memb.comun.para[i]])
	}
			
	ORDER = 1 
	Med.cond.score1 = apply(CondRank1, 2, median)
	Med.cond.score2 = apply(CondRank2, 2, median)

	# CondRank = rbind(c(1:dim(total_cond_scores)[1]), CondRank1, CondRank2)
	# colnames(CondRank) = colnames(data.exp.scale.med)
	# CondRank0 = CondRank[,order(Med.cond.score1, decreasing = T)]

	# CondRank = rbind(c(1:dim(total_cond_scores)[1]), CondRank2, CondRank1)
	# colnames(CondRank) = colnames(data.exp.scale.med)
	# CondRank01 = CondRank[,order(Med.cond.score2, decreasing = T)]

	memb.comun.para.cond = which(rank(Med.cond.score1 - Med.cond.score2, ties.method= "random") < 11)
	memb.comun.cond = which(rank(Med.cond.score1 - Med.cond.score2, ties.method= "random") > 124)
	
	myImagePlot_sym(data.exp.scale[c(memb.comun.para, memb.comun),c(memb.comun.para.cond, memb.comun.cond)])

	aheatmap(data.exp.scale[c(Anch$Ref, Anch$References),c(CondRank0[1,1:10],CondRank01[1,1:10])])
	aheatmap(data.exp.scale[c(Anch$Ref, Anch$References),c(CondRank0[1,1:10],CondRank01[1,1:10])])
	aheatmap(data.exp.scale[Anch$Ref,c(memb.comun.para.cond, memb.comun.cond)])

}



