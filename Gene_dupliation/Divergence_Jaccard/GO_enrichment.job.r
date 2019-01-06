

library(foreach)
library(igraph)
library(csbl.go)
set.prob.table(organism=TAXONOMY.ARABIDOPSIS, type="similarity")
library(doMC)
registerDoMC(30)


## Read condition names
condition.names = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/condition_name.txt", header = T)
NAME = names(table(condition.names))
k = 13
NAME.spe = NAME[k]

## Load generated graph from SA result
load(paste("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Gene_score/Compendium_wise/",NAME[k],"/Graph_V2.RData", sep = ""))

## Read significant pairs
# Result = read.table("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Divergence_Jaccard/Score.txt", header = T)
Result = read.table(paste("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Local/Compendium_wise/",NAME[k],"/Score99/Scores_V2/Anch_Score.txt", sep = ""), header = T)
ref = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/Reference_SD.txt", header = T)
# paralogous = read.table("~/Gene_duplication/Seed_genes/paralogous_total.txt", header = T)
paralogous = read.table("~/Gene_duplication/Seed_genes/paralogous_NoFilter.txt", header = T)
PARA = unique(paralogous$Ref)

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
rownames(RR) = RR$GO


# Add one overlapped genes to all the module's pair in order to consider the cluster size when there is no
# overlapped genes between two modules
Result["Jaccard2"] = (Result$Cluster.overlap + 0.0001)/(Result$Gene.cluster.size + Result$Paralogous.cluster.size)
# Result["D1"] = (Result$Cluster.overlap + 0.0001)/Result$Gene.cluster.size
# Result["D2"] = (Result$Cluster.overlap + 0.0001)/Result$Paralogous.cluster.size
# Result["Dmean"] = apply(Result[, c("D1", "D2")], 1, mean)
# Result["Dsd"] = apply(Result[, c("D1", "D2")], 1, sd)
# Result["Div"] = Result$Dmean/Result$Dsd

##----------------------------------------------------------------------------------------------------------------------------------
## Module GO enrichment analysis
##----------------------------------------------------------------------------------------------------------------------------------

Module.GO = foreach (j = 1:nrow(Result), .combine = rbind) %dopar% {
 
 	print(j)
	
	memb.comun = neighborhood(g, 1, Result$Gene[j], mode = c("out"))[[1]]
	memb.comun.para = neighborhood(g, 1, Result$Paralogous[j], mode = c("out"))[[1]]
	 
	 
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
		
		cat("First modules finished\n")
		names(pvalue) = names(odd) = names(global) = names(loc) = names(freq)
		f = data.frame(pvalue, odd, global, loc)
		f["FDR"] = p.adjust(f$pvalue, method = c("fdr"))
		f = f[order(f$FDR), ]
		top = c(rownames(f[f$FDR < 0.5,]))
		f = f[f$FDR < 0.5,]
		f["GO"] = rownames(f)
		RRf = merge(RR, f, by = "GO")
		RRf = RRf[order(RRf$FDR, decreasing = F), ]
		RRf = cbind.data.frame(RRf, "Jaccard" = Result[j, "Jaccard2"]) 
	} else {
			RRf = data.frame(NA)
		}	
	
	
	## Second module's members by gene name
	members = ref[ref$References %in% memb.comun.para, "Gene"]
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

		cat("Second modules finished\n")
		names(pvalue) = names(odd) = names(global) = names(loc) = names(freq)
		f2 = data.frame(pvalue, odd, global, loc)
		f2["FDR"] = p.adjust(f2$pvalue, method = c("fdr"))
		f2 = f2[order(f2$FDR), ]
		top2 = c(rownames(f2[f2$FDR < 0.5,]))
		f2 = f2[f2$FDR < 0.5,]
		f2["GO"] = rownames(f2)
		RRf2 = merge(RR, f2, by = "GO")
		RRf2 = RRf2[order(RRf2$FDR, decreasing = F), ]
		RRf2 = cbind.data.frame(RRf2, "Jaccard" = Result[j, "Jaccard2"]) 
	} else {
			RRf2 = data.frame(NA)
		}
	
	if(ncol(RRf) > 1 & ncol(RRf2) > 1) rbind(RRf, RRf2)
	if(ncol(RRf) == 1 & ncol(RRf2) > 1) RRf2
	if(ncol(RRf) > 1 & ncol(RRf2) == 1) RRf
	if(ncol(RRf) == 1 & ncol(RRf2) == 1) RRf = RRf
}

write.table(Module.GO, file = paste("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Divergence_Jaccard/development/GO.modules.txt", col.names = T, row.names = F, sep = "\t"))


# --------------------------------------------------------------------------------------------------------------------------------------------------
# For single core jobs, combining all GO enrichment analysis for each module separately
# --------------------------------------------------------------------------------------------------------------------------------------------------
library(preprocessCore)
library(data.table)
setwd("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Divergence_Jaccard/atgenexpress/")
List = list.files()
g = lapply(List, fread)
GO.all = do.call(rbind, g)
# Remove Na values
GO.all = GO.all[!is.na(GO.all$FDR),]

# Only significant GO terms FDR < 0.05
GO.all.sig = GO.all[GO.all$FDR < 0.01,]

Final.score.sig = aggregate(GO.all.sig$Jaccard ~ GO.all.sig$GO, FUN = mean)
Final.score.sig.med = aggregate(GO.all.sig$Jaccard ~ GO.all.sig$GO, FUN = median)
# Final.score.sig.FDR = aggregate(GO.all.sig$FDR ~ GO.all.sig$GO, FUN = median)
Final.score.sig.sd = aggregate(GO.all.sig$Jaccard ~ GO.all.sig$GO, FUN = sd)
AA = table(GO.all.sig$GO)
AA = AA[names(AA) %in% Final.score.sig[,1]]
Final.score.sig = cbind.data.frame(Final.score.sig, "Freq" = AA, "median" = Final.score.sig.med[,2], "SD" = Final.score.sig.sd[,2])

Final.score.sig = Final.score.sig[order(Final.score.sig[, "median"]),]
Final.score.sig["Rank"] = 1:nrow(Final.score.sig)
Final.score.sig.10 = Final.score.sig[Final.score.sig[, "Freq"] > 3, ]
# --------------------------------------------------------------------------------------------------------------------------------------------------



# --------------------------------------------------------------------------------------------------------------------------------------------------
# Permutation test
# --------------------------------------------------------------------------------------------------------------------------------------------------
# For each GO terms and its offsprings (Relationship based on GO trees)
library(GO.db)
library(foreach)
library(doMC)
registerDoMC(3)

# Remove GO terms which are present in the GO.db database
for (i in c(1:nrow(Final.score.sig.10))) {
	print(i)
	get(as.character(Final.score.sig.10[, "GO.all.sig$GO"][i]), GOBPOFFSPRING)
}

Remove.GO = c("GO:0006467", "GO:0006944", "GO:0007243", "GO:0048610")
Final.score.sig.10 = Final.score.sig.10[!(Final.score.sig.10[, "GO.all.sig$GO"] %in% Remove.GO), ]

Group.Ana = foreach (i = c(1:nrow(Final.score.sig.10))) %dopar% {
	print(i)
	group = get(as.character(Final.score.sig.10[, "GO.all.sig$GO"][i]), GOBPOFFSPRING)
	Final.score.sig.10[Final.score.sig.10[, "GO.all.sig$GO"] %in% c(as.character(Final.score.sig.10[, "GO.all.sig$GO"][i]), group),"Rank"]
}

names(Group.Ana) = as.character(Final.score.sig.10[, "GO.all.sig$GO"])
Group.Ana.fin = data.frame("Freq" = sapply(Group.Ana, length), "Average" = sapply(Group.Ana, mean), 
	"Median" = sapply(Group.Ana, median), "SD" = sapply(Group.Ana, sd))
Group.Ana.fin = Group.Ana.fin[Group.Ana.fin$Freq > 1, ]
Group.Ana.fin = Group.Ana.fin[order(Group.Ana.fin$Median),]
Group.Ana.fin["GO"] = rownames(Group.Ana.fin)



GO.size = sort(unique(Group.Ana.fin$Freq))
Permut = matrix(NA, 10000, length(GO.size))
j = 0
for (i in GO.size) {
	print(i)
	j = j + 1
	GO.size.Rand = rep(GO.size[j], 10000)
	# Sample from nodes with the cluster's size
	Permut[, j] = unlist(sapply(GO.size.Rand, function(x){median(Final.score.sig.10[sample(1:nrow(Final.score.sig.10), x), "Rank"])}))
}

Thresh = data.frame("Freq" = GO.size, "Threshold.Div" = apply(Permut, 2, quantile, probs = 0.05), 
	"Threshold.Cons" = apply(Permut, 2, quantile, probs = 0.95))
Group.Ana.fin = merge(Thresh, Group.Ana.fin, by = "Freq")
Group.Ana.fin.Div = Group.Ana.fin[Group.Ana.fin$Median < Group.Ana.fin$Threshold.Div,]
Group.Ana.fin.Cons = Group.Ana.fin[Group.Ana.fin$Median > Group.Ana.fin$Threshold.Cons,]
GO.Div = merge(RR, Group.Ana.fin.Div, by = "GO")
GO.Cons = merge(RR, Group.Ana.fin.Cons, by = "GO")
GO.Div["Pvalue"] = GO.Div$Threshold.Div - GO.Div$Median
GO.Cons["Pvalue"] = GO.Cons$Median - GO.Cons$Threshold.Cons
GO.Div = GO.Div[order(-GO.Div$Pvalue), ]
GO.Cons = GO.Cons[order(-GO.Cons$Pvalue), ]

write.table(GO.Div, file = "/scratch/iwt/ehsab/Gene_duplication/Results/SA/Divergence_Jaccard/GO.Modu.Div.flower.large.txt", quote = F, sep = "\t", row.names = F)
write.table(GO.Cons, file = "/scratch/iwt/ehsab/Gene_duplication/Results/SA/Divergence_Jaccard/GO.Modu.Cons.flower.large.txt", quote = F, sep = "\t", row.names = F)

# --------------------------------------------------------------------------------------------------------------------------------------------------
# Enrichment Map file preparation
# --------------------------------------------------------------------------------------------------------------------------------------------------

# 1. GMT file
library(GO.db)
XX = as.list(GOBPOFFSPRING)

# Remove GO terms which are present in the GO.db database
for (i in c(1:length(XX))) {
	i
	print(get(names(XX[i]), GOBPOFFSPRING))
}

Remove.GO = c("GO:0006467", "GO:0006944", "GO:0007243", "GO:0048610", "GO:0000009")
Final.score.sig.10 = Final.score.sig.10[!(Final.score.sig.10[, "GO.all.sig$GO"] %in% Remove.GO), ]

# Write GO terms which are present in the GO.db database
write.table(NA, file ="/scratch/iwt/ehsab/Gene_duplication/Results/SA/Divergence_Jaccard/abiotic/GMT.txt" , sep = "\t", row.names = F, col.names = F, quote = F, append = F)
for (i in 1:length(XX)) {
	print(i)
	Offspring.GO = get(names(XX[i]), GOBPOFFSPRING)
	if (!is.na(Offspring.GO)) {
		ROW = data.frame(names(XX[i]), RR[names(XX[i]), "X3"], matrix(Offspring.GO, 1, length(Offspring.GO)))
		write.table(ROW, file = "/scratch/iwt/ehsab/Gene_duplication/Results/SA/Divergence_Jaccard/abiotic/GMT.txt", sep = "\t", row.names = F, col.names = F, quote = F, append = T)
	}
}

# 2. Generic result file for GSEA
GO.Div = read.delim("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Divergence_Jaccard/GO.Modu.Div.biotic.large.txt")
GO.Cons = read.delim("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Divergence_Jaccard/GO.Modu.Cons.biotic.large.txt")

GO.Div = merge(RR, GO.Div, by = "GO")
GO.Div["FDR"] = 0.05 - ((GO.Div$Pvalue/GO.Div$Threshold.Div) * 0.0499999999)
GO.Div["P.value"] = GO.Div$FDR - 0.001
GO.Div["Phenotype"] = 1
GO.Cons = merge(RR, GO.Cons, by = "GO")
GO.Cons["FDR"] = 0.05 - ((GO.Cons$Pvalue/GO.Cons$Threshold.Div) * 0.0499999999)
GO.Cons["P.value"] = GO.Cons$FDR - 0.001
GO.Cons["Phenotype"] = -1


Generic.Div = data.frame(GO.Div[,c("GO", "X3.x", "FDR", "P.value", "Phenotype")])
Generic.Cons = data.frame(GO.Cons[,c("GO", "X3.x", "FDR", "P.value", "Phenotype")])
Generic = rbind(Generic.Div, Generic.Cons)
write.table(Generic, file = "/scratch/iwt/ehsab/Gene_duplication/Results/SA/Divergence_Jaccard/GO/EnrichmentMap/flower/Generic.txt", quote = F, row.names = F, sep = "\t")
# --------------------------------------------------------------------------------------------------------------------------------------------------





# --------------------------------------------------------------------------------------------------------------------------------------------------
# Extract significant GO terms from module_GO file and assign Jaccard index to them
# --------------------------------------------------------------------------------------------------------------------------------------------------
library(preprocessCore)
library(data.table)
Score = read.delim("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Combined/Scores.txt")
GOALL = read.delim("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Combined/Module_GO.txt")
GO.all = GOALL[GOALL$Group =="development", ]
# Only significant GO terms FDR < 0.05
GO.all.sig = GO.all[GO.all$FDR < 0.01,]

# Assign Jaccard index to significant GO terms
Score = Score[Score$Group %in% "development",]
Score["Jaccard2"] = (Score$Cluster.overlap + 0.0001)/(Score$Gene.cluster.size + Score$Paralogous.cluster.size)

GO.all.sig = merge(Score, GO.all.sig, by.y = "Seed", by.x = "Gene")

Final.score.sig = aggregate(GO.all.sig$Jaccard2 ~ GO.all.sig$GO, FUN = mean)
Final.score.sig.med = aggregate(GO.all.sig$Jaccard2 ~ GO.all.sig$GO, FUN = median)
# Final.score.sig.FDR = aggregate(GO.all.sig$FDR ~ GO.all.sig$GO, FUN = median)
Final.score.sig.sd = aggregate(GO.all.sig$Jaccard2 ~ GO.all.sig$GO, FUN = sd)
AA = table(GO.all.sig$GO)
AA = AA[names(AA) %in% Final.score.sig[,1]]
Final.score.sig = cbind.data.frame(Final.score.sig, "Freq" = AA, "median" = Final.score.sig.med[,2], "SD" = Final.score.sig.sd[,2])

Final.score.sig = Final.score.sig[order(Final.score.sig[, "median"]),]
Final.score.sig["Rank"] = 1:nrow(Final.score.sig)
Final.score.sig.10 = Final.score.sig[Final.score.sig[, "Freq"] > 3, ]

# Rest of analysis should be followed in the ##For single core jobs, combining all GO enrichment analysis for each module separately## section
# --------------------------------------------------------------------------------------------------------------------------------------------------







