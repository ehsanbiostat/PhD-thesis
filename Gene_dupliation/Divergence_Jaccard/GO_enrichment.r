

library(foreach)
library(igraph)
library(csbl.go)
set.prob.table(organism=TAXONOMY.ARABIDOPSIS, type="similarity")
library(doMC)
registerDoMC(30)


## Read condition names
condition.names = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/condition_name.txt", header = T)
NAME = names(table(condition.names))
k = 2
NAME.spe = NAME[k]

## Load generated graph from SA result
load(paste("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Gene_score/Compendium_wise/",NAME[k],"/Graph_V2.RData", sep = ""))

## Read significant pairs
Result = read.table("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Divergence_Jaccard/Score.txt", header = T)
Result = read.table(paste("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Local/Compendium_wise/",NAME[k],"/Score99/Scores_V2/Anch_Score.txt", sep = ""), header = T)
ref = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/Reference_SD.txt", header = T)
paralogous = read.table("~/Gene_duplication/Seed_genes/paralogous_total.txt", header = T)
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


# Add one overlapped genes to all the module's pair in order to consider the cluster size when there is no
# overlapped genes between two modules
Result["Jaccard2"] = (Result$Cluster.overlap + 0.0001)/(Result$Gene.cluster.size + Result$Paralogous.cluster.size)

# Anchorpoints itself
## ------------------------------------------------------------------------------------------------------------------------------------- ##
Dup.GO = go[go$V3 %in% PARA,]
colnames(Dup.GO)[3] = "Gene"

# Combine anchorpoints and jaccard score into a new data frame
Result.All = data.frame(rbind(data.frame("Gene" = Result$Gene), data.frame("Gene" = Result$Paralogous)), rbind(data.frame(Result$Jaccard2), data.frame(Result$Jaccard2)))
# Select anchorpoints corresponding GO terms
Dup.GO.jac = merge(Result.All, Dup.GO, by = "Gene")

Anch.GO = foreach(i = 1:length(unique(Dup.GO.jac$V2)), .combine = rbind) %do% {
	data.frame("Divergence" = mean(Dup.GO.jac[Dup.GO.jac$V2 %in% unique(Dup.GO.jac$V2)[i], "Result.Jaccard2"]),
				"Divergence.median" = median(Dup.GO.jac[Dup.GO.jac$V2 %in% unique(Dup.GO.jac$V2)[i], "Result.Jaccard2"]),
				"Freq" = length(Dup.GO.jac[Dup.GO.jac$V2 %in% unique(Dup.GO.jac$V2)[i], "Result.Jaccard2"]),
				"SD" = sd(Dup.GO.jac[Dup.GO.jac$V2 %in% unique(Dup.GO.jac$V2)[i], "Result.Jaccard2"]),
				"GO" = unique(Dup.GO.jac$V2)[i])
}


Anch.GO = Anch.GO[order(Anch.GO[, 2], -Anch.GO[, 3]), ]
Anch.GO["Rank"] = 1:nrow(Anch.GO)
write.table(Anch.GO, file = paste("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Divergence_Jaccard/GO.DIV.",NAME.spe,".txt", sep = ""), row.names = F, sep = "\t", quote = F)
## ------------------------------------------------------------------------------------------------------------------------------------- ##



##----------------------------------------------------------------------------------------------------------------------------------
## Permutation test
##----------------------------------------------------------------------------------------------------------------------------------
# For each GO terms and its offsprings (Relationship based on GO trees)
library(GO.db)
library(foreach)
library(doMC)
registerDoMC(3)
# Remove GO terms which are not listed in the GO.db database
Remove.GO = c("GO:0070683", "GO:0006467", "GO:0006944", "GO:0045750", "GO:0007243")
Anch.GO = Anch.GO[!(Anch.GO$GO %in% Remove.GO), ]

for (i in c(1004:nrow(Anch.GO))[-c(128, 589, 725, 943, 1003)]) {
	print(i)
	get(as.character(Anch.GO[, "GO"][i]), GOBPOFFSPRING)
}



Group.Ana = foreach (i = c(1:nrow(Anch.GO))) %dopar% {
	print(i)
	group = get(as.character(Anch.GO[, "GO"][i]), GOBPOFFSPRING)
	Anch.GO[Anch.GO[, "GO"] %in% c(as.character(Anch.GO[, "GO"][i]), group),"Rank"]
}

names(Group.Ana) = as.character(Anch.GO[, "GO"])
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
	Permut[, j] = unlist(sapply(GO.size.Rand, function(x){median(Anch.GO[sample(1:nrow(Anch.GO), x), "Rank"])}))
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

write.table(GO.Div, file = "/scratch/iwt/ehsab/Gene_duplication/Results/SA/Divergence_Jaccard/GO.Anch.Div.abiotic.txt", quote = F, sep = "\t", row.names = F)
write.table(GO.Cons, file = "/scratch/iwt/ehsab/Gene_duplication/Results/SA/Divergence_Jaccard/GO.Anch.Cons.abiotic.txt", quote = F, sep = "\t", row.names = F)



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

