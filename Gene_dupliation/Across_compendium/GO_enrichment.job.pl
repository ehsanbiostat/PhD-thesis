
#!/usr/local/bin/perl -w

use strict;

for(my $i = 1; $i < 4760; $i++){
	
	my $jobfile = "GO".$i.".r";
     open(OUT,">$jobfile");
     print OUT "

library(igraph)
library(foreach)
library(doMC)
# registerDoMC(10)
condition.names = read.table(\"/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/condition_name.txt\", header = T)
NAME = names(table(condition.names))
NAME = NAME[c(1,2,3,6,7,8,9,10,11,12,13,14)]
paralogous = read.table(\"~/Gene_duplication/Seed_genes/paralogous_NoFilter.txt\", header = T)
PARA = unique(paralogous\$Ref)
ref = read.table(\"/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/Reference_SD.txt\", header = T)

# --------------------------------------------------------------------------------------------------------------------------------------------------
# Enriched GO in each module
# --------------------------------------------------------------------------------------------------------------------------------------------------

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

# Module.GO.all = foreach (k = 2:length(NAME), .combine = rbind) %do% {
	k = 1
	load(paste(\"/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Gene_score/Compendium_wise/\",NAME[k],\"/Graph_V2_95.RData\", sep = \"\"))
	# Module.GO = foreach (m = 1:length(PARA), .combine = rbind) %dopar% {
 		m = ".$i."
		memb.comun = neighborhood(g, 1, PARA[m], mode = c(\"out\"))[[1]]
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
			names(pvalue) = names(odd) = names(global) = names(loc) = names(freq)
			f = data.frame(pvalue, odd, global, loc)
			f[\"FDR\"] = p.adjust(f\$pvalue, method = c(\"fdr\"))
			f = f[order(f\$FDR), ]
			top = c(rownames(f[f\$FDR < 0.5,]))
			f = f[f\$FDR < 0.5,]
			f[\"GO\"] = rownames(f)
			RRf = merge(RR, f, by = \"GO\")
			RRf = RRf[order(RRf\$FDR, decreasing = F), ]
			RRf = cbind.data.frame(\"Seed\" = ref[PARA[m], \"Gene\"], RRf, \"Group\" = NAME[k]) 
		} else {
				RRf = as.data.frame(matrix(NA, 1, 12))
				colnames(RRf) = c(\"Seed\", \"GO\", \"X2\", \"X3\", \"X4\", \"X5\", \"pvalue\", \"odd\", \"global\", \"loc\", \"FDR\", \"Group\")
				RRf\$Seed = ref[PARA[m], \"Gene\"]
				RRf\$Group = NAME[k]
			}
		# RRf	
	# }
	write.table(RRf, file = paste(\"/scratch/iwt/ehsab/Gene_duplication/Results/SA/Combined/\",NAME[k],\"\"_95\"/\",m,\".txt\", sep = \"\"), row.names = F, quote = F, sep = \"\\t\")
	# Module.GO
# }

"
	;
	close OUT;
}





