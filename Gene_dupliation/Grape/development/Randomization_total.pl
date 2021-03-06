
#!/usr/local/bin/perl -w

use strict;

for(my $i = 1; $i < 1608; $i++){
		
		my $jobfile = "Random".$i."_99.r";
		 open(OUT,">$jobfile");
		 print OUT "

		 
	
	library(foreach)
	library(doMC)
	library(igraph)
	registerDoMC(8)


	a = read.table(\"/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Compendium_wise/development/Score99/Scores_V2/Score1.txt\", header = T)
	## Select modules having more than 3 members
	cluster.size.Arab = sort(unique(a\$orthologous.cluster.size))[-c(1:3)]

	# nodes
	gene.names = read.table(\"/ngsprojects/iwt/ehsab/Gene_duplication/Integrated CORNET2.0/Vitis_vinifera/gene.names.txt\")
	Ref = read.table(\"/ngsprojects/iwt/ehsab/Gene_duplication/Integrated CORNET2.0/Reference_SD.txt\", header = T)
	load(\"/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Graph_V2.R\")
	map = read.delim(\"/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Orthology/Grape-Arabidopsis/Map.12x.plaza.txt\")

	## Random sampling from anchorpoint nodes
	orthologous = read.table(\"/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Orthology/Grape-Arabidopsis/All.txt\", header = T)

	i = ".$i."
	cluster.size.grape.total = c(4:1607)
		
	cluster.size.grape = rep(cluster.size.grape.total[i], 1000)
	# Sample from nodes with the cluster's size
	cluster.grape = lapply(cluster.size.grape, function(x){sample(as.matrix(gene.names[,1]), x)})
	

	Anch.one = lapply(cluster.grape, function(x) {map[map\$name %in% x, ]})

	Anch.one = lapply(Anch.one, function(x) {x[!(duplicated(x\$name) | duplicated(x\$name, fromLast = TRUE)), \"PLAZA.ID\"]})

	Anch.one = lapply(Anch.one, function(x) {orthologous[orthologous\$PLAZA.ID %in% x, ]})


	Int = function(x, y) {length(unique(x[x\$gene_id.y %in% y, \"PLAZA.ID\"]))}

	Result = foreach (i = 1:length(cluster.size.Arab), .combine=rbind) %dopar%{
		print(i)
		cluster.size.two = rep(cluster.size.Arab[i], 1000)
		Anch.two = lapply(cluster.size.two, function(x){sample(Ref[,1], x)})
		Anc = mapply(Int, Anch.one, Anch.two)
		data.frame(\"Threshold95\" = quantile(Anc, probs = 0.95), \"Threshold99\" = quantile(Anc, probs = 0.99))
		
	}

	write.table(Result, file = paste(\"/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Compendium_wise/development/Score99/Threshold/Threshold\",i,\".txt\", sep = \"\"), col.names = T, row.names = F, sep = \"\\t\")
	"
		;
		close OUT;
	}



