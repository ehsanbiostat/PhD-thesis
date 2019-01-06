
library(igraph)
library(foreach)
library(doMC)
registerDoMC(3)

condition.names = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/condition_name.txt", header = T)
NAME = names(table(condition.names))
k = 1
NAME.spe = NAME[k]
load(paste("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Gene_score/Compendium_wise/",NAME[k],"/Graph_V2.RData", sep = ""))
# paralogous = read.table("~/Gene_duplication/Seed_genes/paralogous_total.txt", header = T)
ref = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/Reference_SD.txt", header = T)
paralogous = read.table("~/Gene_duplication/Seed_genes/paralogous_NoFilter.txt", header = T)
load(paste("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Gene_score/Compendium_wise/",NAME[k],"/Graph_V2.RData", sep = ""))

NO.anch.within = read.delim("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Anch_count/Anch.count.within_1000.txt")
NO.anch.within = NO.anch.within[order(-NO.anch.within$Count), ]
NO.anch.within.10 = NO.anch.within[NO.anch.within$module.size > 9, ]
NO.anch.within.10.20 = NO.anch.within.10[NO.anch.within.10$Count > 20, ]
NO.anch.within.10.20 = NO.anch.within.10.20[order(-NO.anch.within.10.20$Anch.perc), ]
NO.anch.within.10.20.down = NO.anch.within.10[NO.anch.within.10$Count == 0, ]
NO.anch.within.10.20.down = NO.anch.within.10.20.down[order(-NO.anch.within.10.20.down$module.size), ]

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


##----------------------------------------------------------------------------------------------------------------------------------
## Module GO enrichment analysis
##----------------------------------------------------------------------------------------------------------------------------------

Module.GO = foreach (j = 1:20, .combine = rbind) %dopar% {
 
 	print(j)
	
	memb.comun = neighborhood(g, 1, NO.anch.within.10.20.down[j, "PARA"], mode = c("out"))[[1]]
		 
	 
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
		top = c(rownames(f[f$FDR < 0.05,]))
		f = f[f$FDR < 0.05,]
		f["GO"] = rownames(f)
		RRf = merge(RR, f, by = "GO")
		RRf = RRf[order(RRf$FDR, decreasing = F), ]
	} else {
			RRf = data.frame(NA)
		}	
	RRf[1:5,]
}
