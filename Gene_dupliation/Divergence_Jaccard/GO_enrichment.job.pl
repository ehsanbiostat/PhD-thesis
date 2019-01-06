
#!/usr/local/bin/perl -w

use strict;

for(my $i = 1; $i < 5399; $i++){
	
	my $jobfile = "GO".$i.".r";
     open(OUT,">$jobfile");
     print OUT "

library(foreach)
library(igraph)
library(csbl.go)
set.prob.table(organism=TAXONOMY.ARABIDOPSIS, type=\"similarity\")
library(doMC)



## Read condition names
condition.names = read.table(\"/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/condition_name.txt\", header = T)
NAME = names(table(condition.names))
k = 7
NAME.spe = NAME[k]

## Load generated graph from SA result
load(paste(\"/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Gene_score/Compendium_wise/\",NAME[k],\"/Graph_V2.RData\", sep = \"\"))

## Read significant pairs
Result = read.table(paste(\"/scratch/iwt/ehsab/Gene_duplication/Results/SA/Local/Compendium_wise/\",NAME[k],\"/Score99/Scores_V2/Anch_Score.txt\", sep = \"\"), header = T)
ref = read.table(\"/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/Reference_SD.txt\", header = T)
# paralogous = read.table(\"~/Gene_duplication/Seed_genes/paralogous_total.txt\", header = T)
paralogous = read.table(\"~/Gene_duplication/Seed_genes/paralogous_NoFilter.txt\", header = T)
PARA = unique(paralogous\$Ref)

## GO terms loading
go = read.table(\"/ngsprojects/iwt/ehsab/Prioritization/prioritization/Combined_Network/New/GO_ID.txt\",sep = \"\\t\")
clas = c(\"C\", \"F\", \"P\")
k = 3
go = go[go[, 4] == clas[k], ]

goo = readLines(\"/ngsprojects/iwt/ehsab/Prioritization/prioritization/Combined_Network/TAIR/GO_terms.txt\")
RR = matrix(NA, 1, 5)
for (i in 1:length(goo)){
	RR = rbind(RR, strsplit(goo[i], \"\\t\")[[1]])
}
RR = RR[-1,]
RR = data.frame(RR)
colnames(RR)[1] = \"GO\"


# Add one overlapped genes to all the module's pair in order to consider the cluster size when there is no
# overlapped genes between two modules
Result[\"Jaccard2\"] = (Result\$Cluster.overlap + 0.0001)/(Result\$Gene.cluster.size + Result\$Paralogous.cluster.size)


##----------------------------------------------------------------------------------------------------------------------------------
## Module GO enrichment analysis
##----------------------------------------------------------------------------------------------------------------------------------

# Module.GO = foreach (j = 1:nrow(Result), .combine = rbind) %dopar% {
 
 	j = ".$i."
	
	memb.comun = neighborhood(g, 1, Result\$Gene[j], mode = c(\"out\"))[[1]]
	memb.comun.para = neighborhood(g, 1, Result\$Paralogous[j], mode = c(\"out\"))[[1]]
	 
	 
	## First module's members by gene name
	members = ref[ref\$References %in% memb.comun, \"Gene\"]
	 ## Extract GO terms
	go.terms = go[go\$V1 %in% members, \"V2\"]

	freq = table(matrix(go.terms)) # Frequency of GO among
	pvalue = odd = global = loc = c()
	
	if (length(go.terms) > 0){
		for (i in 1:length(freq)){
			z = matrix(c(length(go.terms) - freq[i], freq[i]), 2)
			x = matrix(table(go[,2] == names(freq[i])))
			x = cbind(x ,z) 
			pvalue[i] = fisher.test(x, alternative = \"greater\")\$p
			global[i] = x[2,1]
			loc[i] = x[2,2]
			odd[i] = (x[1,1] * x[2,2])/(x[1,2] * x[2,1])
		}
		
		cat(\"First modules finished\n\")
		names(pvalue) = names(odd) = names(global) = names(loc) = names(freq)
		f = data.frame(pvalue, odd, global, loc)
		f[\"FDR\"] = p.adjust(f\$pvalue, method = c(\"fdr\"))
		f = f[order(f\$FDR), ]
		top = c(rownames(f[f\$FDR < 0.5,]))
		f = f[f\$FDR < 0.5,]
		f[\"GO\"] = rownames(f)
		RRf = merge(RR, f, by = \"GO\")
		RRf = RRf[order(RRf\$FDR, decreasing = F), ]
		RRf = cbind.data.frame(RRf, \"Jaccard\" = Result[j, \"Jaccard2\"]) 
	} else {
			RRf = data.frame(NA)
		}	
	
	
	## Second module's members by gene name
	members = ref[ref\$References %in% memb.comun.para, \"Gene\"]
	go.terms = go[go\$V1 %in% members, \"V2\"]

	freq = table(matrix(go.terms)) # Frequency of GO among
	pvalue = odd = global = loc = c()

	if (length(go.terms) > 0){
		for (i in 1:length(freq)){
			z = matrix(c(length(go.terms) - freq[i], freq[i]), 2)
			x = matrix(table(go[,2] == names(freq[i])))
			x = cbind(x ,z) 
			pvalue[i] = fisher.test(x, alternative = \"greater\")\$p
			global[i] = x[2,1]
			loc[i] = x[2,2]
			odd[i] = (x[1,1] * x[2,2])/(x[1,2] * x[2,1])
		}

		cat(\"Second modules finished\n\")
		names(pvalue) = names(odd) = names(global) = names(loc) = names(freq)
		f2 = data.frame(pvalue, odd, global, loc)
		f2[\"FDR\"] = p.adjust(f2\$pvalue, method = c(\"fdr\"))
		f2 = f2[order(f2\$FDR), ]
		top2 = c(rownames(f2[f2\$FDR < 0.5,]))
		f2 = f2[f2\$FDR < 0.5,]
		f2[\"GO\"] = rownames(f2)
		RRf2 = merge(RR, f2, by = \"GO\")
		RRf2 = RRf2[order(RRf2\$FDR, decreasing = F), ]
		RRf2 = cbind.data.frame(RRf2, \"Jaccard\" = Result[j, \"Jaccard2\"]) 
	} else {
			RRf2 = data.frame(NA)
		}
	
	if(ncol(RRf) > 1 & ncol(RRf2) > 1) Module.GO = rbind(RRf, RRf2)
	if(ncol(RRf) == 1 & ncol(RRf2) > 1) Module.GO = RRf2
	if(ncol(RRf) > 1 & ncol(RRf2) == 1) Module.GO = RRf
	if(ncol(RRf) == 1 & ncol(RRf2) == 1) Module.GO = RRf = RRf
# }

write.table(Module.GO, file = paste(\"/scratch/iwt/ehsab/Gene_duplication/Results/SA/Divergence_Jaccard/\",NAME.spe,\"/GO.modules\",j,\".txt\", sep = \"\"), col.names = T, row.names = F, sep = \"\\t\", quote = F)

"
	;
	close OUT;
}





