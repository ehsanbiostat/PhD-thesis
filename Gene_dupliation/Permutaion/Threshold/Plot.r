
paralogous = read.table("~/Gene_duplication/Seed_genes/paralogous_total.txt", header = T)
p = sort(unique(c(as.matrix(paralogous))))


setwd("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Gene_score/Permutation/SA_threshold/")
list = read.table("~/Gene_duplication/Scripts/General/order.SA_threshold.txt")
list = as.vector(unlist(list))

g = lapply(list, function(x) mget(load(x)))


AA = function(x){length(which(gene_score_module[, x] > 0))}
Len1 = sapply(p, AA)
load("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Gene_score/Permutation/SA_threshold/Graph1.6.RData")
Len2 = sapply(p, AA)

load("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Gene_score/Permutation/SA_threshold/Graph1.7.RData")
Len3 = sapply(p, AA)

load("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Gene_score/Permutation/SA_threshold/Graph1.8.RData")
Len4 = sapply(p, AA)

load("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Gene_score/Permutation/SA_threshold/Graph1.9.RData")
Len5 = sapply(p, AA)

load("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Gene_score/Permutation/SA_threshold/Graph2.RData")
Len6 = sapply(p, AA)

load("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Gene_score/Permutation/SA_threshold/Graph2.1.RData")
Len7 = sapply(p, AA)

load("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Gene_score/Permutation/SA_threshold/Graph2.2.RData")
Len8 = sapply(p, AA)

load("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Gene_score/Permutation/SA_threshold/Graph2.3.RData")
Len9 = sapply(p, AA)

load("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Gene_score/Permutation/SA_threshold/Graph2.4.RData")
Len10 = sapply(p, AA)

Len = list(Len1, Len2, Len3, Len4, Len5, Len6, Len7, Len8, Len9, Len10)
med = sapply(Len, median)/median(Len1)

jpeg("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Gene_score/Permutation/SA_threshold/Changing.jpeg", width = 12, height = 9, units = "in", pointsize = 12, res = 500)
plot(sapply(Len, median))
dev.off()

















