
##----------------------------------------------------------------------------------------------------------------------------------
## Shared anchorpoints GO enrichment analysis
##----------------------------------------------------------------------------------------------------------------------------------

library(foreach)
library(igraph)
library(csbl.go)
set.prob.table(organism=TAXONOMY.ARABIDOPSIS, type="similarity")
library(doMC)
registerDoMC(1)


## Read condition names
condition.names = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated CORNET2.0/condition_name.txt", header = T)
NAME = names(table(condition.names))
NAME.spe = NAME[6]

## Load generated graph from SA result
load(paste("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Gene_score/Compendium_wise/",NAME.spe,"/Graph_V2.RData", sep = ""))

## Read significant pairs
ref = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated CORNET2.0/Reference_SD.txt", header = T)



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



GO.BP = GO.BP.Anch = list()

# SemSim.Observe = foreach (j = 1:nrow(Result), .combine = c) %dopar% {
	j = "4567_15903"
	cut.off = 0.05
	 print(j)
	 Anchor.shared = read.table(paste("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Local/Compendium_wise/",NAME.spe,"/Score99/Module_Anchorpoints_V2/Module",j,".txt", sep = ""), header = T)	
	 members = Anchor.shared$Gene.y
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
		 top = c(rownames(f[f$FDR < cut.off,]))
		 f = f[f$FDR < cut.off,]
		 f["GO"] = rownames(f)
		 RRf = merge(RR, f, by = "GO")
		 RRf = RRf[order(RRf$FDR, decreasing = F), ]
	} else {
			RRf1 = data.frame(NA)
			top1 = "NA"
			}	
	
	 members = Anchor.shared$Gene.x
	 go.terms = go[go$V1 %in% members, "V2"]
	if (length(go.terms) > 0){
		 freq = table(matrix(go.terms)) # Frequency of GO among
		 pvalue = odd = global = loc = c()

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
		 top2 = c(rownames(f2[f2$FDR < cut.off, ]))
		 f2 = f2[f2$FDR < cut.off,]
		 f2["GO"] = rownames(f2)
		 RRf2 = merge(RR, f2, by = "GO")
		 RRf2 = RRf2[order(RRf2$FDR, decreasing = F), ]
	} else {
			RRf2 = data.frame(NA)
			top2 = "NA"
			}
	 
	 write.table(RRf, file = paste("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Local/Validation/Compendium_wise/",NAME.spe,"/V2/GO/Anchorpoints/Enriched_terms/GO.terms",j,"_1.txt", sep = "") 
	 , append = F, col.names = T, row.names = F, sep = "\t", quote = F)
	 
	 write.table(RRf2, file = paste("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Local/Validation/Compendium_wise/",NAME.spe,"/V2/GO/Anchorpoints/Enriched_terms/GO.terms",j,"_2.txt", sep = "") 
	 , append = F, col.names = T, row.names = F, sep = "\t", quote = F)
	 
	 
# }



write.table(SemSim.Observe, file = paste("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Local/Validation/Compendium_wise/",NAME.spe,"/V2/GO/Anchorpoints/GO.Total.txt", sep = ""), 
col.names = F, row.names = F, sep = "\t")
 
semsim = c()
SemSim.Randomization = foreach (j = 1:2000, .combine = c) %dopar% {
	# print(j)
	a = readLines(paste("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Local/Validation/Compendium_wise/",NAME.spe,"/V2/GO/Anchorpoints/",j,".txt", sep = ""))
	N = length(strsplit(a[1], "\t")[[1]]) - 1
	M = length(strsplit(a[2], "\t")[[1]]) - 1
	
	for (i in 1:10){
		test1 = t(as.matrix(c(as.matrix("A"),as.vector(sample(go$V2, N, replace = F)))))
		test2 = t(as.matrix(c(as.matrix("B"),as.vector(sample(go$V2, M, replace = F)))))
			
		write.table(test1, file = paste("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Local/Validation/Compendium_wise/",NAME.spe,"/V2/GO/test/",j,".txt", sep="") 
		, append = F, col.names = F, row.names = F, sep = "\t", quote = F)

		write.table(test2, file = paste("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Local/Validation/Compendium_wise/",NAME.spe,"/V2/GO/test/",j,".txt", sep="") 
		, append = T, col.names = F, row.names = F, sep = "\t", quote = F)
			
		ent = entities.from.text(paste("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Local/Validation/Compendium_wise/",NAME.spe,"/V2/GO/test/",j,".txt", sep=""))
		semsim[i] = entity.sim.many(ent, "BP", "ResnikGraSM")[1,2]
	}
	
	quantile(semsim, probs = 0.95, na.rm = T)
	
}

write.table(SemSim.Randomization, file = paste("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Local/Validation/Compendium_wise/",NAME.spe,"/V2/GO/Anchorpoints/SemSim.Randomization.txt", sep = ""), col.names = F , row.names = F, sep = "\t")



 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
