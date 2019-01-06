
library(foreach)
# install.packages("doMC", repos = getOption("http://cran.us.r-project.org"))
library(doMC)
registerDoMC(cores = 10)
# install.packages("doMPI", repos = getOption("http://cran.us.r-project.org"))
# library(doMPI)
# cl = startMPIcluster(16)
# registerDoMPI(cl)

load("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Compendium_wise/development/Score99/Threshold/Threshold.RData")


condition.names = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated CORNET2.0/condition_name.txt", header = T)
NAME = names(table(condition.names))


k = 6

Result = foreach (i = 1:19589) %dopar%{
	print(i)
	a = read.table(paste("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Compendium_wise/",NAME[k],"/Score99/Scores_V2/Score",i,".txt", sep = ""), header = T)
	# Select clusters with more than 3 nodes
	a1 = a[a$orthologous.cluster.size > 3,]
	
	# Select relevant threshold (Normalized_Score & Jaccard) values for corresponding Gene and paralogouse's cluster size
	bb = Threshold[Threshold$Grape == a1$Gene.cluster.size[1], ]
	colnames(bb)[2] = c("orthologous.cluster.size")
	

	a11 = merge(bb, a1, by = c("orthologous.cluster.size"))
	# Select those pairs which have higher paralogous's score value than threshold & less Jaccard's score value than threshold
	Result = a11[(a11$orthologous.cluster > a11$Threshold95) & a11$orthologous.cluster > 1, ]
	# inter = intersect(a12$Paralogous, paralogous[paralogous$Ref %in% i, "Anch"])
	# Average = mean(a12[a12$Paralogous %in% inter, "Paralogous.cluster"])
	# Result[i] = length(intersect(a12$Paralogous, paralogous[paralogous$Ref %in% i, "Anch"]))
	# data.frame(nrow(a12), print(Result[i]), print(length(paralogous[paralogous$Ref %in% i, "Anch"])), Average)
	write.table(Result, file = paste("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Compendium_wise/",NAME[k],"/Score99/Result95_V2_",i,".txt", sep = ""), row.names = F, col.names = T, sep = "\t")
}

# write.table(Result, file = paste("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Compendium_wise/",NAME[k],"/Score99/Result95_V2.txt", sep = ""), row.names = F, col.names = T, sep = "\t")


















