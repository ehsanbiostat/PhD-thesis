
library(igraph)
library(foreach)
library(doMC)
registerDoMC(5)
condition.names = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/condition_name.txt", header = T)
NAME = names(table(condition.names))
NAME = NAME[c(1,2,3,6,7,8,9,10,11,12,13,14)]
paralogous = read.table("~/Gene_duplication/Seed_genes/paralogous_NoFilter.txt", header = T)
PARA = unique(paralogous$Ref)
ref = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/Reference_SD.txt", header = T)

# --------------------------------------------------------------------------------------------------------------------------------------------------
# Scores
# --------------------------------------------------------------------------------------------------------------------------------------------------
Result = foreach (k = 1:length(NAME), .combine = rbind) %dopar% {
	print(k)
	load(paste("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Gene_score/Compendium_wise/",NAME[k],"/Graph_V2.RData", sep = ""))
	Score = matrix(NA, nrow(paralogous), 9)
	colnames(Score) = c("Gene", "Paralogous", "Paralogous.cluster", "Gene.cluster.size", "Paralogous.cluster.size", "Cluster.overlap", 
		"Paralogous.withi", "WGDevent", "KS")
	Score = as.data.frame(Score)
	for (m in 1:nrow(paralogous)) {
		# print(m)
		i = paralogous[m,1]
		j = paralogous[m,2]
		memb.comun = data.frame("Ref" = neighborhood(g, 1, i, mode = c("out"))[[1]])
		memb.comun.para = merge(memb.comun, paralogous, by = "Ref")
		memb.comun.para.size = nrow(memb.comun.para)
		Score[m, "Gene"] = i
		Score[m, "Paralogous"] = j
		Score[m, "Paralogous.cluster"] = length(intersect(neighborhood(g, 1, j, mode = c("out"))[[1]], memb.comun.para$Anch))
		Score[m, "Paralogous.withi"] = nrow(paralogous[paralogous$Ref %in% neighborhood(g, 1, j, mode = c("out"))[[1]],])
		Score[m, "Cluster.overlap"] = length(intersect(neighborhood(g, 1, j, mode = c("out"))[[1]], neighborhood(g, 1, i, mode = c("out"))[[1]]))
		Score[m, "Gene.cluster.size"] = dim(memb.comun)[1]
		Score[m, "Paralogous.cluster.size"] = length(neighborhood(g, 1, j, mode = c("out"))[[1]])
		Score[m, "WGDevent"] = as.vector(paralogous[m, "WGDevent"])
		Score[m, "KS"] = paralogous[m,"KS"]
		}
	Score["Normalized.Paralogous.cluster"] = Score$Paralogous.cluster/(Score$Gene.cluster.size + Score$Paralogous.cluster.size)
	Score["Normalized.Paralogous.cluster.Gene"] = Score$Paralogous.cluster/memb.comun.para.size
	Score["Normalized.Paralogous.whitin.Gene"] = memb.comun.para.size/Score$Gene.cluster.size
	Score["Normalized.Paralogous.cluster.Paralogous"] = Score$Paralogous.cluster/Score$Paralogous.withi
	Score["Normalized.Paralogous.whitin.Paralogous"] = Score$Paralogous.withi/Score$Paralogous.cluster.size
	Score["Jaccard"] = Score$Cluster.overlap/(Score$Gene.cluster.size + Score$Paralogous.cluster.size)
	Score["Jaccard2"] = (Score$Cluster.overlap + 0.0001)/(Score$Gene.cluster.size + Score$Paralogous.cluster.size)
	Score[is.na(Score)] = 0
	write.table(Score, file = paste("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Local/Compendium_wise/",NAME[k],"/Score99/Scores_V2/Anch_Score.txt", sep = ""), row.names = F, quote = F, sep = "\t")
	Score["Group"] = NAME[k]
	# Score = Score[order(Score$Jaccard2),]
	Score["Rank"] = 1:nrow(Score)
	Score
}

Result$Gene = ref[Result$Gene, "Gene"]
Result$Paralogous = ref[Result$Paralogous, "Gene"]
write.table(Result, file = "/scratch/iwt/ehsab/Gene_duplication/Results/SA/Combined/Scores.txt", row.names = F, quote = F, sep = "\t")
# --------------------------------------------------------------------------------------------------------------------------------------------------




# --------------------------------------------------------------------------------------------------------------------------------------------------
# Membership
# --------------------------------------------------------------------------------------------------------------------------------------------------
Membership = foreach (k = 1:length(NAME), .combine = rbind) %dopar% {
	print(k)
	load(paste("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Gene_score/Compendium_wise/",NAME[k],"/Graph_V2.RData", sep = ""))
	Memb = foreach (m = 1:length(PARA), .combine = rbind) %do% {
		vicinity = neighborhood(g, 1, PARA[m], mode = c("out"))[[1]]
		data.frame("Seed" = vicinity[1], "Members" = vicinity)[-1,]
	}
	data.frame(Memb, "Group" = NAME[k])
}


Membership$Seed = ref[Membership$Seed, "Gene"]
Membership$Members = ref[Membership$Members, "Gene"]
write.table(Membership, file = "/scratch/iwt/ehsab/Gene_duplication/Results/SA/Combined/Membership.txt", row.names = F, quote = F, sep = "\t")
# --------------------------------------------------------------------------------------------------------------------------------------------------




# --------------------------------------------------------------------------------------------------------------------------------------------------
# Enriched GO in each module
# --------------------------------------------------------------------------------------------------------------------------------------------------

## GO terms loading
go = read.table("/ngsprojects/iwt/ehsab/Prioritization/prioritization/Combined_Network/New/GO_ID.txt",sep = "\t")
clas = c("C", "F", "P")
k = 3
go = go[go[, 4] == clas[k], ]

goo = readLines("/ngsprojects/iwt/ehsab/Prioritization/prioritization/Combined_Network/TAIR/GO_terms.txt")
RR = matrix(NA, 1, 5)
for (i in 1:length(goo)){
	RR = rbind(RR, strsplit(goo[i], "\t")[[1]])
}
RR = RR[-1,]
RR = data.frame(RR)
colnames(RR)[1] = "GO"

Module.GO.all = foreach (k = 1:length(NAME), .combine = rbind) %do% {
	load(paste("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Gene_score/Compendium_wise/",NAME[k],"/Graph_V2.RData", sep = ""))
	Module.GO = foreach (j = 1:length(PARA), .combine = rbind) %dopar% {
 		print(j)
		memb.comun = neighborhood(g, 1, PARA[m], mode = c("out"))[[1]]
		## First module's members by gene name
		members = ref[ref$References %in% memb.comun, "Gene"]
		 ## Extract GO terms
		go.terms = go[go$V1 %in% members, "V2"]
		freq = table(matrix(go.terms)) # Frequency of GO among
		pvalue = odd = global = loc = c()
		if (length(go.terms) > 0){
			for (i in 1:length(freq)){
				z = matrix(c(length(go.terms) - freq[i], freq[i]), 2)
				x = matrix(table(go[,2] == names(freq[i])))
				x = cbind(x ,z) 
				pvalue[i] = fisher.test(x, alternative = "greater")$p
				global[i] = x[2,1]
				loc[i] = x[2,2]
				odd[i] = (x[1,1] * x[2,2])/(x[1,2] * x[2,1])
			}
			names(pvalue) = names(odd) = names(global) = names(loc) = names(freq)
			f = data.frame(pvalue, odd, global, loc)
			f["FDR"] = p.adjust(f$pvalue, method = c("fdr"))
			f = f[order(f$FDR), ]
			top = c(rownames(f[f$FDR < 0.5,]))
			f = f[f$FDR < 0.5,]
			f["GO"] = rownames(f)
			RRf = merge(RR, f, by = "GO")
			RRf = RRf[order(RRf$FDR, decreasing = F), ]
			RRf = cbind.data.frame("Seed" = ref[PARA[m], "Gene"], RRf, "Group" = NAME[k]) 
		} else {
				RRf = data.frame(NA)
			}
		RRf	
	}
	write.table(Module.GO, file = paste("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Combined/Module_GO_",NAME[k],".txt", sep = ""), row.names = F, quote = F, sep = "\t")
	Module.GO
}

write.table(Module.GO.all, file = "/scratch/iwt/ehsab/Gene_duplication/Results/SA/Combined/Module_GO.txt", col.names = T, row.names = F, sep = "\t"))
# --------------------------------------------------------------------------------------------------------------------------------------------------




# --------------------------------------------------------------------------------------------------------------------------------------------------
# CORNET
# --------------------------------------------------------------------------------------------------------------------------------------------------


Module.GO.all = foreach (k = 1:length(NAME), .combine = rbind) %do% {
# Condition scores
load(paste("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Condition_score/Compendium_wise/",NAME[k],"/Raw/C",PARA[m],".RData", sep = ""))

# Top 10 conditions for up-regulated & down-regulated
up = names(cond_scores[order(-cond_scores),])
down = names(cond_scores[order(cond_scores),])

# top up-regulated conditions
# Remove "X" and "." from condition names
conditions = as.numeric(unlist(strsplit(sub("X", "", up), "[.]")))
conditions = conditions[conditions > 4]
# Add seed gene name, directionality, compendium and PO annotaion of condition for top up-regulated conditions
UP = data.frame("Seed" = ref[PARA[m], "Gene"], t(data.frame(COR.annot.PO[names(COR.annot.PO) %in% conditions])), 
	"CORNT" = conditions, "Direction" = "Up", "Rank" = 1:length(up), "Group" = NAME[k])

# top down-regulated conditions
# Remove "X" and "." from condition names
conditions = as.numeric(unlist(strsplit(sub("X", "", down), "[.]")))
conditions = conditions[conditions > 4]
# Add seed gene name, directionality, compendium and PO annotaion of condition
DOWN = data.frame("Seed" = ref[PARA[m], "Gene"], t(data.frame(COR.annot.PO[names(COR.annot.PO) %in% conditions])), 
	"CORNT" = conditions, "Direction" = "Down", "Rank" = 1:length(up), "Group" = NAME[k])

# Concatenate UP and DOWN conditions
rbind(UP, DOWN)



# --------------------------------------------------------------------------------------------------------------------------------------------------
Agg = cbind(aggregate(Result$Rank ~ Result$Gene + Result$Paralogous, FUN = median), aggregate(Result$Rank ~ Result$Gene + Result$Paralogous, FUN = sd), 
	aggregate(Result$Cluster.overlap ~ Result$Gene + Result$Paralogous, FUN = function(x){sum(x == 0)}))
Agg = Agg[, -c(4:5, 7:8)]
colnames(Agg)[1] = "Ref"
Agg = merge(Agg, paralogous, by = "Ref")
Agg = Agg[, -2]
colnames(Agg)[2:4] = c("median", "SD", "Count")
Agg = Agg[order(-Agg$Count, Agg$median, Agg$SD),]
head(Agg)
Result[Result$Gene %in% Agg[1, "Ref"], c("Cluster.overlap","Gene.cluster.size", "Paralogous.cluster.size", "Group", "Rank")]

Com.wise = list()
Com.wise = foreach (k = 1:length(NAME), .combine = rbind) %dopar% {
for (k in 1:length(NAME)) {
	print(k)
	load(paste("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Gene_score/Compendium_wise/",NAME[k],"/Graph_V2.RData", sep = ""))
	go.terms = go[go$V3 %in% unlist(neighborhood(g, 1, Agg[1, "Ref"], mode = c("out"))), "V2"]

	# GO enrichment analysis
	freq = table(matrix(go.terms)) # Frequency of GO among
	pvalue = odd = global = loc = c()
	if (length(go.terms) > 0){
		for (i in 1:length(freq)){
			z = matrix(c(length(go.terms) - freq[i], freq[i]), 2)
			x = matrix(table(go[,2] == names(freq[i])))
			x = cbind(x ,z) 
			pvalue[i] = fisher.test(x, alternative = "greater")$p
			global[i] = x[2,1]
			loc[i] = x[2,2]
			odd[i] = (x[1,1] * x[2,2])/(x[1,2] * x[2,1])
		}
		names(pvalue) = names(odd) = names(global) = names(loc) = names(freq)
		f = data.frame(pvalue, odd, global, loc)
		f["FDR"] = p.adjust(f$pvalue, method = c("fdr"))
		f = f[order(f$FDR), ]
		top = c(rownames(f[f$FDR < 0.9,]))
		f = f[f$FDR < 0.9,]
		f["GO"] = rownames(f)
		RRf = merge(RR, f, by = "GO")
		RRf = RRf[order(RRf$FDR, decreasing = F), ]
		RRf = cbind.data.frame(RRf, "Group" = NAME[k], "Ref" = Agg[1, "Ref"])
	} else {
			RRf = data.frame(NA)
		}	
	

	go.terms = go[go$V3 %in% unlist(neighborhood(g, 1, Agg[1, "Anch"], mode = c("out"))), "V2"]
	freq = table(matrix(go.terms)) # Frequency of GO among
	pvalue = odd = global = loc = c()

	if (length(go.terms) > 0){
		for (i in 1:length(freq)){
			z = matrix(c(length(go.terms) - freq[i], freq[i]), 2)
			x = matrix(table(go[,2] == names(freq[i])))
			x = cbind(x ,z) 
			pvalue[i] = fisher.test(x, alternative = "greater")$p
			global[i] = x[2,1]
			loc[i] = x[2,2]
			odd[i] = (x[1,1] * x[2,2])/(x[1,2] * x[2,1])
		}

		names(pvalue) = names(odd) = names(global) = names(loc) = names(freq)
		f2 = data.frame(pvalue, odd, global, loc)
		f2["FDR"] = p.adjust(f2$pvalue, method = c("fdr"))
		f2 = f2[order(f2$FDR), ]
		top2 = c(rownames(f2[f2$FDR < 0.9,]))
		f2 = f2[f2$FDR < 0.9,]
		f2["GO"] = rownames(f2)
		RRf2 = merge(RR, f2, by = "GO")
		RRf2 = RRf2[order(RRf2$FDR, decreasing = F), ]
		RRf2 = cbind.data.frame(RRf2, "Group" = NAME[k], "Ref" = Agg[1, "Anch"])
	} else {
			RRf2 = data.frame(NA)
		}
	
	if(ncol(RRf) > 1 & ncol(RRf2) > 1) Com.wise[k] = list(rbind(RRf, RRf2))
	if(ncol(RRf) == 1 & ncol(RRf2) > 1) Com.wise[k] = list(RRf2)
	if(ncol(RRf) > 1 & ncol(RRf2) == 1) Com.wise[k] = list(RRf)
	# if(ncol(RRf) == 1 & ncol(RRf2) == 1) RRf = RRf
}



# --------------------------------------------------------------------------------------------------------------------------------------------------
# Module size distribution
# --------------------------------------------------------------------------------------------------------------------------------------------------
library(data.table)
library(beanplot)
Module.size = foreach (k = 1:length(NAME), .combine = rbind) %dopar% {
	print(k)
	load(paste("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Gene_score/Compendium_wise/",NAME[k],"/Graph_V2.RData", sep = ""))
	data.frame("Size" = degree(g, mode = "out"), "Group" = NAME[k])
}	

Module.size.col = foreach (k = 1:length(NAME), .combine = cbind) %dopar% {
	print(k)
	load(paste("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Gene_score/Compendium_wise/",NAME[k],"/Graph_V2.RData", sep = ""))
	degree(g, mode = "out")
}	



Module.out = foreach (k = 1:length(NAME)) %dopar% {
	print(k)
	load(paste("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Gene_score/Compendium_wise/",NAME[k],"/Graph_V2.RData", sep = ""))
	neighborhood(g, 1, PARA, mode = "out")
}	

# Intersection for each seed gene across all 12 compendia
AA = matrix(NA, length(PARA), length(NAME) ^ 2)
for(k in 1:length(PARA)) {
	print(k)
	BB = sapply(Module.out, function(x){x[[k]]})
	Jac = foreach (i = 1:12, .combine = c) %do% {
		mapply(function(x,y){length(intersect(x,y))/length(union(x,y))}, BB, BB[i])
	}
	AA[k, ] = Jac
}

# Module size distribution
jpeg("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Across_compendia/Plot/ModuleSizeBeanPlot.jpeg", width=8, height=5, units="in", res=100, pointsize = 5)
ggplot(Module.size, aes(factor(Group), Size)) +  geom_violin(fill = "grey80", colour = "#3366FF", adjust = 0.5, trim = FALSE, scale = "width") + geom_boxplot(width=.1)
dev.off()

Result = fread("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Combined/Scores.txt")
GO = fread("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Combined/Module_GO.txt")
Memb = fread("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Combined/Membership.txt")
# --------------------------------------------------------------------------------------------------------------------------------------------------



sort(table(Memb[Memb[, Seed == ref[PARA[1],1]], Members]), decreasing= F)




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



#################################################################################################################################################################################
# Gene score correlation for two example across all 14 compendia for selected anchorpoints vs other genes within module
#################################################################################################################################################################################
library(ggplot2)
library(foreach)

condition.names = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/condition_name.txt", header = T)
NAME = names(table(condition.names))[-c(4,5)]

Anchorpoint = read.delim("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/Anchorpoint/Geneset1_2_53_90.txt.txt")
Anchorpoint1 = read.delim("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/Anchorpoint/AT1G50090_AT3G19710_52_90.txt")

ref = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/Reference_SD.txt", header = T)
rownames(ref) = ref[,1]

Leaf = read.table("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/Module/Leaf_90.txt")
Pollen = read.table("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/Module/Pollen_90.txt")

PollenRef = ref[as.character(Pollen[,1]), 2]
LeafRef = ref[as.character(Leaf[,1]), 2]

PollenAnch = Anch$Pollen
PollenAnch = Anch$Root
LeafAnch = Anch$Leaf

PollenLeaf = foreach(j = 1:12, .combine = cbind) %do% {
	load(paste("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Gene_score/Compendium_wise/",as.name(NAME[j]),"/Raw/G",PollenRef[i],".RData", sep = ""))
	rbind(merge(a, Anchorpoint, by.x = c("Var1", "Var2"), by.y = c("Pollen", "Leaf")), merge(a, Anchorpoint, by.x = c("Var2", "Var1"), by.y = c("Pollen", "Leaf")))[, "Freq"]
}

RootLeaf = foreach(j = 1:12, .combine = cbind) %do% {
	a = read.delim(paste("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Gene_score/Compendium_wise/",as.name(NAME[j]),"/Gene.score.correlation.txt", sep = ""))
	rbind(merge(a, Anchorpoint1, by.x = c("Var1", "Var2"), by.y = c("Root", "Leaf")), merge(a, Anchorpoint1, by.x = c("Var2", "Var1"), by.y = c("Root", "Leaf")))[, "Freq"]
}

colnames(PollenLeaf) = colnames(RootLeaf) = NAME

jpeg(paste("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Gene_score/Compendium_wise/development/PollenLeaf_GeneScoreCorrelation.jpeg", sep = ""), width = 11, height = 7, units = "in", res = 500, pointsize=8)
myImagePlot(PollenLeaf[,c(names(sort(apply(PollenLeaf, 2, median))))])
dev.off()

jpeg(paste("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Gene_score/Compendium_wise/development/RootLeaf_GeneScoreCorrelation.jpeg", sep = ""), width = 11, height = 7, units = "in", res = 500, pointsize=8)
myImagePlot(RootLeaf[,c(names(sort(apply(RootLeaf, 2, median))))])
dev.off()


RootLeafTotal = foreach(j = 1:12, .combine = cbind) %do% {
	read.delim(paste("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Gene_score/Compendium_wise/",as.name(NAME[j]),"/Gene.score.correlation.txt", sep = ""))[, "Freq"]
}

colnames(PollenLeaf) = colnames(RootLeafTotal) = NAME

#################################################################################################################################################################################












#################################################################################################################################################################################
# Gene score correlation for two example across all 14 compendia for selected anchorpoints
#################################################################################################################################################################################
library(ggplot2)
library(foreach)

condition.names = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/condition_name.txt", header = T)
NAME = names(table(condition.names))[-c(4,5)]

Anchorpoint = read.delim("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/Anchorpoint/Geneset1_2_53_90.txt.txt")
Anchorpoint1 = read.delim("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/Anchorpoint/AT1G50090_AT3G19710_52_90.txt")


PollenLeaf = foreach(j = 1:12, .combine = cbind) %do% {
	a = read.delim(paste("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Gene_score/Compendium_wise/",as.name(NAME[j]),"/Gene.score.correlation.txt", sep = ""))
	rbind(merge(a, Anchorpoint, by.x = c("Var1", "Var2"), by.y = c("Pollen", "Leaf")), merge(a, Anchorpoint, by.x = c("Var2", "Var1"), by.y = c("Pollen", "Leaf")))[, "Freq"]
}

RootLeaf = foreach(j = 1:12, .combine = cbind) %do% {
	a = read.delim(paste("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Gene_score/Compendium_wise/",as.name(NAME[j]),"/Gene.score.correlation.txt", sep = ""))
	rbind(merge(a, Anchorpoint1, by.x = c("Var1", "Var2"), by.y = c("Root", "Leaf")), merge(a, Anchorpoint1, by.x = c("Var2", "Var1"), by.y = c("Root", "Leaf")))[, "Freq"]
}

colnames(PollenLeaf) = colnames(RootLeaf) = NAME

jpeg(paste("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Gene_score/Compendium_wise/development/PollenLeaf_GeneScoreCorrelation.jpeg", sep = ""), width = 11, height = 7, units = "in", res = 500, pointsize=8)
myImagePlot(PollenLeaf[,c(names(sort(apply(PollenLeaf, 2, median))))])
dev.off()

jpeg(paste("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Gene_score/Compendium_wise/development/RootLeaf_GeneScoreCorrelation.jpeg", sep = ""), width = 11, height = 7, units = "in", res = 500, pointsize=8)
myImagePlot(RootLeaf[,c(names(sort(apply(RootLeaf, 2, median))))])
dev.off()


RootLeafTotal = foreach(j = 1:12, .combine = cbind) %do% {
	read.delim(paste("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Gene_score/Compendium_wise/",as.name(NAME[j]),"/Gene.score.correlation.txt", sep = ""))[, "Freq"]
}

colnames(PollenLeaf) = colnames(RootLeafTotal) = NAME

#################################################################################################################################################################################










#################################################################################################################################################################################
# Gene score correlation distribution for diverged APs compare to the all APs
#################################################################################################################################################################################
library(ggplot2)
library(foreach)
options(width = 200)
condition.names = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/condition_name.txt", header = T)
NAME = names(table(condition.names))[-c(4,5)]
source("/ngsprojects/iwt/ehsab/Gene_duplication/Scripts/General/Multiple_plot.r")

Anchorpoint = read.delim("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/Anchorpoint/Geneset1_2_53_90.txt.txt")
Anchorpoint1 = read.delim("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/Anchorpoint/AT1G50090_AT3G19710_52_90.txt")


for(j in 1:12) { # Four graphs in one page, density comparison as well and violion plot for gene score correlation distribution
	print(j)
	a = read.delim(paste("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Gene_score/Compendium_wise/",as.name(NAME[j]),"/Gene.score.correlation.txt", sep = ""))
	b = rbind(merge(a, Anchorpoint, by.x = c("Var1", "Var2"), by.y = c("Pollen", "Leaf")), merge(a, Anchorpoint, by.x = c("Var2", "Var1"), by.y = c("Pollen", "Leaf")))
	a["Group"] = "All anchorpoints"
	aa = a[(a$References.x %in% b$References.x & a$References.y %in% b$References.y),]
	aa[, "Group"] = "Pollen-Leaf"
	a = rbind(a, aa)
	
	b = rbind(merge(a, Anchorpoint1, by.x = c("Var1", "Var2"), by.y = c("Root", "Leaf")), merge(a, Anchorpoint1, by.x = c("Var2", "Var1"), by.y = c("Root", "Leaf")))
	aa = a[(a$References.x %in% b$References.x & a$References.y %in% b$References.y),]
	aa[, "Group"] = "Root-Leaf"
	a = rbind(a, aa)
	

	aAncient = a[a$WGDevent == "ancient" & !(is.na(a$WGDevent)),]
	aRecent = a[a$WGDevent == "recent" & !(is.na(a$WGDevent)),]

		
	AnPo = round(wilcox.test(a[a$WGDevent == "ancient" & !(is.na(a$WGDevent)) & a$Group == "All anchorpoints", "Freq"], a[a$WGDevent == "ancient" & !(is.na(a$WGDevent)) & a$Group == "Pollen-Leaf", "Freq"], alternative = c("greater"))$p.value, 4)
	AnRo = round(wilcox.test(a[a$WGDevent == "ancient" & !(is.na(a$WGDevent)) & a$Group == "All anchorpoints", "Freq"], a[a$WGDevent == "ancient" & !(is.na(a$WGDevent)) & a$Group == "Root-Leaf", "Freq"], alternative = c("greater"))$p.value, 4)
	

	RePo = round(wilcox.test(a[a$WGDevent == "recent" & !(is.na(a$WGDevent)) & a$Group == "All anchorpoints", "Freq"], a[a$WGDevent == "recent" & !(is.na(a$WGDevent)) & a$Group == "Pollen-Leaf", "Freq"], alternative = c("greater"))$p.value, 4)
	ReRo =round(wilcox.test(a[a$WGDevent == "recent" & !(is.na(a$WGDevent)) & a$Group == "All anchorpoints", "Freq"], a[a$WGDevent == "recent" & !(is.na(a$WGDevent)) & a$Group == "Root-Leaf", "Freq"], alternative = c("greater"))$p.value, 4)
	
			
	p1 = ggplot(aAncient, aes(x = Freq, fill = Group)) +
	 geom_density(alpha = 0.3) + ggtitle("Ancient event") + 
	 scale_fill_discrete(name="Distribution_P-value", breaks=c("All anchorpoints", "Pollen-Leaf", "Root-Leaf"), labels=c("All anchorpoints", paste("Pollen-Leaf", AnPo, sep = "_"), paste("Root-Leaf", AnRo, sep = "_"))) + xlab("Gene score correlation") + theme(text = element_text(size=8))

	p2 = ggplot(aRecent, aes(x = Freq, fill = Group)) +
	 geom_density(alpha = 0.3) + ggtitle("Recent event") + 
	 scale_fill_discrete(name="Distribution_P-value", breaks=c("All anchorpoints", "Pollen-Leaf", "Root-Leaf"), labels=c("All anchorpoints", paste("Pollen-Leaf", RePo, sep = "_"), paste("Root-Leaf", ReRo, sep = "_"))) + xlab("Gene score correlation") + theme(text = element_text(size=8))

	P3 = ggplot(aAncient, aes(y = Freq, x = factor(Group))) +
	 geom_violin(fill = "lightgreen") +
	 geom_boxplot(width=.2, color = "black", fill = "yellow") +
	 ylab("Gene score correlation") + xlab("") + ggtitle("Ancient event") + theme(text = element_text(size=8)) + 
	 annotate("text", x = "Pollen-Leaf", y = -1, label = paste("Pvalue=",AnPo, sep = ""), size = 2.5) + annotate("text", x = "Root-Leaf", y = -1, label = paste("Pvalue=",AnRo, sep = ""), size = 2.5)
	
	P4 = ggplot(aRecent, aes(y = Freq, x = factor(Group))) +
	 geom_violin(fill = "lightgreen") +
	 geom_boxplot(width=.2, color = "black", fill = "yellow") +
	 ylab("Gene score correlation") + xlab("") + ggtitle("Recent event") + theme(text = element_text(size=8)) +
	 annotate("text", x = "Pollen-Leaf", y = -1, label = paste("Pvalue=",RePo, sep = ""), size = 2.5) + annotate("text", x = "Root-Leaf", y = -1, label = paste("Pvalue=",ReRo, sep = ""), size = 2.5)


	jpeg(paste("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Gene_score/Across_compendia/GeneScoreCorrDist/",as.name(NAME[j]),".jpeg", sep = ""), width = 15, height = 8, units = "in", res = 800, pointsize=8)
	multiplot(p1, p2,P3, P4, cols = 2)
	dev.off()

}

#################################################################################################################################################################################



#################################################################################################################################################################################
# Gene score correlation overlap for stress, hormone and development 
#################################################################################################################################################################################

library(ggplot2)
library(foreach)
options(width = 200)
source("/ngsprojects/iwt/ehsab/Gene_duplication/Scripts/General/Multiple_plot.r")
condition.names = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/condition_name.txt", header = T)
NAME = names(table(condition.names))[-c(4,5)]

Anchorpoint = read.delim("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/Anchorpoint/Geneset1_2_53_90.txt.txt")
Anchorpoint1 = read.delim("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/Anchorpoint/AT1G50090_AT3G19710_52_90.txt")
GeneScoreCorr = read.delim("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Gene_score/Across_compendia/GeneScoreCorr.txt")

GeneScoreCorrAncient = GeneScoreCorr[GeneScoreCorr$WGDevent == "ancient" & !(is.na(GeneScoreCorr$WGDevent)),]
GeneScoreCorrRecent = GeneScoreCorr[GeneScoreCorr$WGDevent == "recent" & !(is.na(GeneScoreCorr$WGDevent)),]


# Pollen = rbind(merge(GeneScoreCorr, Anchorpoint, by.x = c("Var1", "Var2"), by.y = c("Pollen", "Leaf")), merge(GeneScoreCorr, Anchorpoint, by.x = c("Var2", "Var1"), by.y = c("Pollen", "Leaf")))
Pollen = GeneScoreCorr[(GeneScoreCorr$Var1 %in% Anchorpoint$Pollen & GeneScoreCorr$Var2 %in% Anchorpoint$Leaf) | (GeneScoreCorr$Var1 %in% Anchorpoint$Leaf & GeneScoreCorr$Var2 %in% Anchorpoint$Pollen), ]
Root = GeneScoreCorr[(GeneScoreCorr$Var1 %in% Anchorpoint1$Root & GeneScoreCorr$Var2 %in% Anchorpoint1$Leaf) | (GeneScoreCorr$Var1 %in% Anchorpoint1$Leaf & GeneScoreCorr$Var2 %in% Anchorpoint1$Root), ]

PollenAncient = Pollen[Pollen$WGDevent == "ancient" & !(is.na(Pollen$WGDevent)),]
PollenRecent = Pollen[Pollen$WGDevent == "recent" & !(is.na(Pollen$WGDevent)),]
RootAncient = Root[Root$WGDevent == "ancient" & !(is.na(Root$WGDevent)),]
RootRecent = Root[Root$WGDevent == "recent" & !(is.na(Root$WGDevent)),]

scale_colour_manual("Homeologs", values = c("Pollen" = "red", "Root" = "darkgreen")) +
	  scale_linetype_manual(name="Legend", values = c("Pollen" = "dotted", "Root" = "dotted")) +
	  theme(axis.text = element_text(color = "black"),legend.key.height  = grid::unit(0.1, "npc"))


for(i in 4:15) {
	p1 = ggplot(GeneScoreCorr, aes(x = GeneScoreCorr[, i])) + geom_histogram(fill = "white", color = "blue", binwidth = 0.01) +
	 geom_vline(xintercept = Pollen[,i], color = "red", show_guide = TRUE) + geom_vline(xintercept = Root[,i], color = "darkgreen", show_guide = TRUE) + xlab(NAME[i-3]) + ggtitle("Gene score correlation for homeologs") +
	  theme(text = element_text(size=8))

	p2 = ggplot(GeneScoreCorrAncient, aes(x = GeneScoreCorrAncient[,i])) + geom_histogram(fill = "white", color = "blue", binwidth = 0.01) +
	 geom_vline(xintercept = PollenAncient[,i], color = "red") + geom_vline(xintercept = RootAncient[,i], color = "darkgreen") + xlab("Beta WGDevent") +
	  theme(text = element_text(size=8))

	p3 = ggplot(GeneScoreCorrRecent, aes(x = GeneScoreCorrRecent[,i])) + geom_histogram(fill = "white", color = "blue", binwidth = 0.01) +
	 geom_vline(xintercept = PollenRecent[,i], color = "red") + geom_vline(xintercept = RootRecent[,i], color = "darkgreen") + xlab("Alpha WGDevent") +
	  theme(text = element_text(size=8))

	jpeg(paste("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Gene_score/Across_compendia/",NAME[i-3],".jpeg", sep = ""), width = 11, height = 7, units = "in", res = 500, pointsize=8)
	multiplot(p1, p2, p3,  cols = 1)
	dev.off()
}
#################################################################################################################################################################################


gplot <- ggplot(df, aes(x=1:length(error), y=error)) +
  scale_x_continuous(breaks = seq_along(df$error)) +
  geom_point(size = 3) +
  geom_line() +
  geom_errorbar(data = df, aes(x = 1:length(error), ymin = error - sDev, ymax = error + sDev), 
                width = 0.1) + 
  geom_hline(data = df, aes(yintercept = error[minimum] + sDev[minimum], linetype="a", colour="a")) +
  geom_vline(data= data.frame(type="b", col="b", minimum=minimum), 
             aes(linetype=type, colour=col, xintercept = minimum), size = 1, show_guide = TRUE) +
  geom_vline(data= data.frame(type="b", col="b", best.model=best.model),
             aes(linetype="c", colour="c", xintercept = best.model), size = 1, show_guide = TRUE) +
  scale_colour_manual(name="Legend", values = c("a" = "black", "b" = "red", "c" = "blue")) +
  scale_linetype_manual(name="Legend", values = c("a" = "dashed", "b" = "dotted", "c" = "dotted")) +
  theme_gray(base_size = 18) + 
  theme(axis.text = element_text(color = "black"),
        legend.key.height  = grid::unit(0.1, "npc")) +
  labs(x = "# of parameters", fontface = "bold") + 
  labs(y = "CV error") +
  labs(title = "Cross-validation error curve")



i = 4
p4 = ggplot(GeneScoreCorr, aes(x = GeneScoreCorr[, 4])) + geom_histogram(fill = "white", color = "blue", binwidth = 0.01) +
	  geom_vline(xintercept = Root[,i], color = "darkgreen") + ggtitle(NAME[i-3]) +
	  theme(text = element_text(size=8)) + coord_flip() + scale_y_reverse() + xlab("")
i = 5
p5 = ggplot(GeneScoreCorr, aes(x = GeneScoreCorr[, 5])) + geom_histogram(fill = "white", color = "blue", binwidth = 0.01) +
	  geom_vline(xintercept = Root[,i], color = "darkgreen") + ggtitle(NAME[i-3]) +
	  theme(text = element_text(size=8)) + coord_flip() + scale_y_reverse() + xlab("")
i = 6
p6 = ggplot(GeneScoreCorr, aes(x = GeneScoreCorr[, 6])) + geom_histogram(fill = "white", color = "blue", binwidth = 0.01) +
	  geom_vline(xintercept = Root[,i], color = "darkgreen") + ggtitle(NAME[i-3]) +
	  theme(text = element_text(size=8)) + coord_flip() + scale_y_reverse() + xlab("")
i = 7
p7 = ggplot(GeneScoreCorr, aes(x = GeneScoreCorr[, 7])) + geom_histogram(fill = "white", color = "blue", binwidth = 0.01) +
	  geom_vline(xintercept = Root[,i], color = "darkgreen") + ggtitle(NAME[i-3]) +
	  theme(text = element_text(size=8)) + coord_flip() + scale_y_reverse() + xlab("")
i = 8
p8 = ggplot(GeneScoreCorr, aes(x = GeneScoreCorr[, 8])) + geom_histogram(fill = "white", color = "blue", binwidth = 0.01) +
	  geom_vline(xintercept = Root[,i], color = "darkgreen") + ggtitle(NAME[i-3]) +
	  theme(text = element_text(size=8)) + coord_flip() + scale_y_reverse() + xlab("")
i = 9
p9 = ggplot(GeneScoreCorr, aes(x = GeneScoreCorr[, 9])) + geom_histogram(fill = "white", color = "blue", binwidth = 0.01) +
	  geom_vline(xintercept = Root[,i], color = "darkgreen") + ggtitle(NAME[i-3]) +
	  theme(text = element_text(size=8)) + coord_flip() + scale_y_reverse() + xlab("")
i = 10
p10 = ggplot(GeneScoreCorr, aes(x = GeneScoreCorr[, 10])) + geom_histogram(fill = "white", color = "blue", binwidth = 0.01) +
	  geom_vline(xintercept = Root[,i], color = "darkgreen") + ggtitle(NAME[i-3]) +
	  theme(text = element_text(size=8)) + coord_flip() + scale_y_reverse() + xlab("")
i = 11
p11 = ggplot(GeneScoreCorr, aes(x = GeneScoreCorr[, 11])) + geom_histogram(fill = "white", color = "blue", binwidth = 0.01) +
	  geom_vline(xintercept = Root[,i], color = "darkgreen") + ggtitle(NAME[i-3]) +
	  theme(text = element_text(size=8)) + coord_flip() + scale_y_reverse() + xlab("")
i = 12
p12 = ggplot(GeneScoreCorr, aes(x = GeneScoreCorr[, 12])) + geom_histogram(fill = "white", color = "blue", binwidth = 0.01) +
	  geom_vline(xintercept = Root[,i], color = "darkgreen") + ggtitle(NAME[i-3]) +
	  theme(text = element_text(size=8)) + coord_flip() + scale_y_reverse() + xlab("")
i = 13
p13 = ggplot(GeneScoreCorr, aes(x = GeneScoreCorr[, 13])) + geom_histogram(fill = "white", color = "blue", binwidth = 0.01) +
	  geom_vline(xintercept = Root[,i], color = "darkgreen") + ggtitle(NAME[i-3]) +
	  theme(text = element_text(size=8)) + coord_flip() + scale_y_reverse() + xlab("")
i = 14
p14 = ggplot(GeneScoreCorr, aes(x = GeneScoreCorr[, 14])) + geom_histogram(fill = "white", color = "blue", binwidth = 0.01) +
	  geom_vline(xintercept = Root[,i], color = "darkgreen") + ggtitle(NAME[i-3]) +
	  theme(text = element_text(size=8)) + coord_flip() + scale_y_reverse() + xlab("")
i = 15
p15 = ggplot(GeneScoreCorr, aes(x = GeneScoreCorr[, 15])) + geom_histogram(fill = "white", color = "blue", binwidth = 0.01) +
	  geom_vline(xintercept = Root[,i], color = "darkgreen") + ggtitle(NAME[i-3]) +
	  theme(text = element_text(size=8)) + coord_flip() + scale_y_reverse() + xlab("")

jpeg(paste("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Gene_score/Across_compendia/Across_RootLeaf.jpeg", sep = ""), width = 15, height = 8, units = "in", res = 500, pointsize=8)
multiplot(p7, p4, p14, p6, p10, p15, p12, p13, p11, p8, p5, p9,  cols = 12)
dev.off()
apply(Root[, 4:15], 2, function(x) {sum(x < -0.5)})
apply(Root[, 4:15], 2, median)



i = 4
p4 = ggplot(GeneScoreCorr, aes(x = GeneScoreCorr[, 4])) + geom_histogram(fill = "white", color = "blue", binwidth = 0.01) +
	  geom_vline(xintercept = Pollen[,i], color = "darkgreen") + ggtitle(NAME[i-3]) +
	  theme(text = element_text(size=8)) + coord_flip() + scale_y_reverse() + xlab("")
i = 5
p5 = ggplot(GeneScoreCorr, aes(x = GeneScoreCorr[, 5])) + geom_histogram(fill = "white", color = "blue", binwidth = 0.01) +
	  geom_vline(xintercept = Pollen[,i], color = "darkgreen") + ggtitle(NAME[i-3]) +
	  theme(text = element_text(size=8)) + coord_flip() + scale_y_reverse() + xlab("")
i = 6
p6 = ggplot(GeneScoreCorr, aes(x = GeneScoreCorr[, 6])) + geom_histogram(fill = "white", color = "blue", binwidth = 0.01) +
	  geom_vline(xintercept = Pollen[,i], color = "darkgreen") + ggtitle(NAME[i-3]) +
	  theme(text = element_text(size=8)) + coord_flip() + scale_y_reverse() + xlab("")
i = 7
p7 = ggplot(GeneScoreCorr, aes(x = GeneScoreCorr[, 7])) + geom_histogram(fill = "white", color = "blue", binwidth = 0.01) +
	  geom_vline(xintercept = Pollen[,i], color = "darkgreen") + ggtitle(NAME[i-3]) +
	  theme(text = element_text(size=8)) + coord_flip() + scale_y_reverse() + xlab("")
i = 8
p8 = ggplot(GeneScoreCorr, aes(x = GeneScoreCorr[, 8])) + geom_histogram(fill = "white", color = "blue", binwidth = 0.01) +
	  geom_vline(xintercept = Pollen[,i], color = "darkgreen") + ggtitle(NAME[i-3]) +
	  theme(text = element_text(size=8)) + coord_flip() + scale_y_reverse() + xlab("")
i = 9
p9 = ggplot(GeneScoreCorr, aes(x = GeneScoreCorr[, 9])) + geom_histogram(fill = "white", color = "blue", binwidth = 0.01) +
	  geom_vline(xintercept = Pollen[,i], color = "darkgreen") + ggtitle(NAME[i-3]) +
	  theme(text = element_text(size=8)) + coord_flip() + scale_y_reverse() + xlab("")
i = 10
p10 = ggplot(GeneScoreCorr, aes(x = GeneScoreCorr[, 10])) + geom_histogram(fill = "white", color = "blue", binwidth = 0.01) +
	  geom_vline(xintercept = Pollen[,i], color = "darkgreen") + ggtitle(NAME[i-3]) +
	  theme(text = element_text(size=8)) + coord_flip() + scale_y_reverse() + xlab("")
i = 11
p11 = ggplot(GeneScoreCorr, aes(x = GeneScoreCorr[, 11])) + geom_histogram(fill = "white", color = "blue", binwidth = 0.01) +
	  geom_vline(xintercept = Pollen[,i], color = "darkgreen") + ggtitle(NAME[i-3]) +
	  theme(text = element_text(size=8)) + coord_flip() + scale_y_reverse() + xlab("")
i = 12
p12 = ggplot(GeneScoreCorr, aes(x = GeneScoreCorr[, 12])) + geom_histogram(fill = "white", color = "blue", binwidth = 0.01) +
	  geom_vline(xintercept = Pollen[,i], color = "darkgreen") + ggtitle(NAME[i-3]) +
	  theme(text = element_text(size=8)) + coord_flip() + scale_y_reverse() + xlab("")
i = 13
p13 = ggplot(GeneScoreCorr, aes(x = GeneScoreCorr[, 13])) + geom_histogram(fill = "white", color = "blue", binwidth = 0.01) +
	  geom_vline(xintercept = Pollen[,i], color = "darkgreen") + ggtitle(NAME[i-3]) +
	  theme(text = element_text(size=8)) + coord_flip() + scale_y_reverse() + xlab("")
i = 14
p14 = ggplot(GeneScoreCorr, aes(x = GeneScoreCorr[, 14])) + geom_histogram(fill = "white", color = "blue", binwidth = 0.01) +
	  geom_vline(xintercept = Pollen[,i], color = "darkgreen") + ggtitle(NAME[i-3]) +
	  theme(text = element_text(size=8)) + coord_flip() + scale_y_reverse() + xlab("")
i = 15
p15 = ggplot(GeneScoreCorr, aes(x = GeneScoreCorr[, 15])) + geom_histogram(fill = "white", color = "blue", binwidth = 0.01) +
	  geom_vline(xintercept = Pollen[,i], color = "darkgreen") + ggtitle(NAME[i-3]) +
	  theme(text = element_text(size=8)) + coord_flip() + scale_y_reverse() + xlab("") + scale_colour_manual("", breaks = c("TempMax", "TempMedia", "TempMin"),values = c("red", "green", "blue"))

jpeg(paste("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Gene_score/Across_compendia/Across_PollenLeaf.jpeg", sep = ""), width = 15, height = 8, units = "in", res = 500, pointsize=8)
multiplot(p7, p4, p12, p10, p14, p5, p13, p15, p11, p9, p6, p8,  cols = 12)
dev.off()

apply(Pollen[, 4:15], 2, function(x) {sum(x < -0.5)})
apply(Pollen[, 4:15], 2, median)




#################################################################################################################################################################################
# Gene score correlation overlap for stress, hormone and development 
#################################################################################################################################################################################
library(ggplot2)
library(foreach)

condition.names = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/condition_name.txt", header = T)
NAME = names(table(condition.names))[-c(4,5)]

Anchorpoint = read.delim("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/Anchorpoint/Geneset1_2_53_90.txt.txt")
Anchorpoint1 = read.delim("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/Anchorpoint/AT1G50090_AT3G19710_52_90.txt")

a = read.delim(paste("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Gene_score/Compendium_wise/",as.name(NAME[j]),"/Gene.score.correlation.txt", sep = ""))
abiotic = a[order(a$Freq),]

GeneScoreCorr = foreach(j = 1:12, .combine = cbind) %do% {
	a = read.delim(paste("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Gene_score/Compendium_wise/",as.name(NAME[j]),"/Gene.score.correlation.txt", sep = ""))
	AA = data.frame(a$Freq)
	colnames(AA) = NAME[j]
	AA
}

GeneScoreCorr = cbind(a[,1:3], GeneScoreCorr)

rbind(merge(GeneScoreCorr, Anchorpoint1, by.x = c("Var1", "Var2"), by.y = c("Root", "Leaf")), merge(GeneScoreCorr, Anchorpoint1, by.x = c("Var2", "Var1"), by.y = c("Root", "Leaf")))
rbind(merge(GeneScoreCorr, Anchorpoint, by.x = c("Var1", "Var2"), by.y = c("Pollen", "Leaf")), merge(GeneScoreCorr, Anchorpoint, by.x = c("Var2", "Var1"), by.y = c("Pollen", "Leaf")))[, -c(16:19)]


AA = as.character(GeneScoreCorr[order(GeneScoreCorr$abiotic), "name.combine"][1:500])
BB = as.character(GeneScoreCorr[order(GeneScoreCorr$development), "name.combine"][1:500])
CC = as.character(GeneScoreCorr[order(GeneScoreCorr$hormone), "name.combine"][1:500])

colnames(PollenLeaf) = colnames(RootLeaf) = NAME

jpeg(paste("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Gene_score/Compendium_wise/development/PollenLeaf_GeneScoreCorrelation.jpeg", sep = ""), width = 11, height = 7, units = "in", res = 500, pointsize=8)
myImagePlot(PollenLeaf[,c(names(sort(apply(PollenLeaf, 2, median))))])
dev.off()

jpeg(paste("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Gene_score/Compendium_wise/development/RootLeaf_GeneScoreCorrelation.jpeg", sep = ""), width = 11, height = 7, units = "in", res = 500, pointsize=8)
myImagePlot(RootLeaf[,c(names(sort(apply(RootLeaf, 2, median))))])
dev.off()

#################################################################################################################################################################################




