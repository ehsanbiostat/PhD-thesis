library(foreach)
library(doMC)
library(preprocessCore)
registerDoMC(16)

load("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Compendium_wise/abiotic/Score99/Result/Total1.RData")

# ICC results
load("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Compendium_wise/abiotic/ICC.RData")

# Select paralogous modules and their ancesteral modules from grapvine
# A pair of paralogous modules should have a same grapevine modules which has a significant enrichment score
A = sapply(g, unlist)
AA = sapply(A, length)

# Load paralogous modules from Arabidopsis data
load("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Compendium_wise/abiotic/Score99/Result.RData")

# Load significant grape modules modules with arabidopsis modules
load("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Compendium_wise/abiotic/Score99/G.score.RData")

load("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Compendium_wise/abiotic/Score99/Result/Intersect_Orth.RData")

# Repeat row function
rep.row <- function(x,n){
   matrix(rep(x,each=n),nrow=n)
}


## Grape
	library(igraph)
	condition.names = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated CORNET2.0/condition_name.txt", header = T)
	NAME = names(table(condition.names))

	k = 1
	
	## Grape modules
	load("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Compendium_wise/stress/Graph_V2.R")
	g.grape = g
	
	## Arabidopsis modules
	load(paste("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Gene_score/Compendium_wise/",NAME[k],"/Graph_V2.RData", sep = ""))
	
	orthologous = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Orthology/Grape-Arabidopsis/All.txt", header = T)
	gene.names = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated CORNET2.0/Vitis_vinifera/Stress/gene.names.txt")
	Ref = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated CORNET2.0/Reference_SD.txt", header = T)
	paralogous = read.table("~/Gene_duplication/Seed_genes/paralogous_total.txt", header = T)


# Quantile normalization of all orthology scores
TT.scale = normalize.quantiles(G.score[, -c(1:9)])
colnames(TT.scale) = colnames(G.score)[-c(1:9)]
G.score.copy = G.score
G.score = cbind(TT.scale, "orthologous" = G.score[, "orthologous"], 
	"Gene" = G.score[, "Gene"])


# Quantile normalization of all paralogy scores
TT.scale.Para = normalize.quantiles(Result[, -c(1:12)])
colnames(TT.scale.Para) = colnames(Result)[-c(1:12)]
Result.copy = Result
Result = cbind(TT.scale.Para, "Paralogous" = Result[, "Paralogous"], 
	"Gene" = Result[, "Gene"])




TT = foreach (m = 1:length(which(AA > 1)), .combine = rbind) %dopar% { 
	print(m)
	
	# Extract all potential pre-duplicated modules for each paralogous module
	i = unlist(A[which(AA > 1)[m]])
	# Extract Arabidopsis modules which shared potential pre-duplicated modules
	j = Result[which(AA > 1)[m], "Gene"]
	k = Result[which(AA > 1)[m], "Paralogous"]


	# Select orthologous scores
	graph.line = data.frame(G.score[G.score[, "Gene"] %in% i & (G.score[, "orthologous"] == j | G.score[, "orthologous"] == k), c("Normalized.orthologous.cluster", "Normalized.orthologous.cluster.grape" ,"Normalized.orthologous.cluster.Arab" , "Normalized.orthologous.whitin.grape" , "Normalized.orthologous.whitin.Arab" , "orthologous", "Gene")])

	# Select paralogous scores
	Arab.scores = data.frame(rep.row(Result[which(AA > 1)[m], c(1:6)], nrow(graph.line)/2))
	colnames(Arab.scores) = colnames(Result)[1:6]

	# ICC score
	CORR = foreach (p = 1:length(Correlation_total[[m]]), .combine = rbind) %do% {
		sapply(Correlation_total[[m]][[p]], median)
	}
	colnames(CORR) = c("Cor1", "Cor2")

	# Generate row for each Arabidopsis-Grapevine module combination and all
	# related scores
	cbind(graph.line[seq.int(1, nrow(graph.line), 2), ], graph.line[seq.int(2, nrow(graph.line), 2), ], data.frame(INTER.Tot1[m]), data.frame(CORR), Arab.scores)
}


save(TT, file = "/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Compendium_wise/abiotic/Combined_Score.RData")

# TT = TT[!duplicated(TT[,c(6,13,14)]),]
# TT = TT[-which(is.na(TT[, 18])), ]
# TT = TT[-which(is.na(TT[, 19])), ]
TT["Cor.Diff"] = abs(TT$Cor1 - TT$Cor2)
TTscore = apply(TT[, -c(6:7,13:19, 25:26)], 1, sum)
TT["Score"] = TTscore
# TT = TT[order(-TT$Score),]
TT.no.zero = TT[TT$INTER > 0, ]
TT.no.zero = TT.no.zero[order(-TT.no.zero$INTER), ]
TT.filter = TT.no.zero[which(TT.no.zero$INTER == 5), ]
TT.filter["Cor.sum"] = TT.filter$Cor1 + TT.filter$Cor2
TT.filter["COR"] = TT.filter$Cor.Diff + TT.filter$Cor.sum
TT.filter = TT.filter[order(TT.filter[, "Score"] ,TT.filter[, "COR"],  decreasing = T), ]
TT.para = TT.filter[, c(6, 13)]
TT.para = TT.para[!duplicated(TT.para), ]
write.table(TT.para, file = "/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Validation/Compendium_wise/development/1.Para.list.txt", row.names = F, quote = F, sep = "\t")

TT.para.no = TT[, c("Gene", "Gene.1")]
TT.para.no = TT.para.no[!duplicated(TT.para.no), ]
write.table(TT.para.no, file = "/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Validation/Compendium_wise/abiotic/Para.list.txt", row.names = F, quote = F, sep = "\t")


# Extract set genes
TT.para.orth = TT.filter[, c(6, 13:14)]
TT.para.orth = TT.para.orth[!duplicated(TT.para.orth), ]
write.table(TT.para.orth, file = "/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Validation/Compendium_wise/development/5.Para.Orth.list.txt", row.names = F, quote = F, sep = "\t")





TT.no.zero = TT.no.zero[order(TT.no.zero[, "INTER"] ,TT.no.zero[, "Cor.Diff"],  decreasing = T), ]

List = c("orthologous.1", "Gene.1", "orthologous", "Score", "Cor1", "Cor2", "Cor.Diff")
fit.zero = prcomp(TT.filter[, List[-c(1:3)]], retx = TRUE, center = T, scale = T)
varimax3.zero <- varimax(fit.zero$rotation[, 1:2])
New.score = fit.zero$x %*% varimax3.zero$loadings

TT.filter = TT.filter[order(TT.filter[, "Score"] ,TT.filter[, "Cor.Diff"],  decreasing = T), ]

fit.zero = prcomp(TT.no.zero[, -c(6:7, 13:14, 16:17)], retx = TRUE)
varimax3.zero <- varimax(fit.zero$rotation[,1:3])


# library(rgl)
# plot3d(fit$scores[,1:3])


library(scatterplot3d)
SCORE = data.frame(newData <- as.matrix(TT.filter[, List[-c(1:3)]]) %*% varimax3.zero$loadings)

with(SCORE, {
   s3d = scatterplot3d(PC1, PC2, PC3,        # x y and z axis
                 color="blue", pch=19, # filled blue circles
                 type="h",             # lines to the horizontal plane
                 main="3-D Scatterplot Example 2",
                 xlab="PC1",
                 ylab="PC2",
                 zlab="PC3")
   s3d.coords <- s3d$xyz.convert(PC1, PC2, PC3) # convert 3D coords to 2D projection
    text(s3d.coords$x, s3d.coords$y,             # x and y coordinates
         labels=row.names(SCORE),               # text to plot
         cex=.5, pos=4)
})


with(TT.filter, {
   s3d = scatterplot3d(Cor1, Score, Cor2,        # x y and z axis
                 color="blue", pch=19, # filled blue circles
                 type="h",             # lines to the horizontal plane
                 main="3-D Scatterplot Example 2",
                 xlab="COR1",
                 ylab="Score",
                 zlab="COR2")
   s3d.coords <- s3d$xyz.convert(Cor1, Score, Cor2) # convert 3D coords to 2D projection
    text(s3d.coords$x, s3d.coords$y,             # x and y coordinates
         labels=row.names(TT.filter),               # text to plot
         cex=.5, pos=4)
})


> plot(TT.filter[, "Score"], TT.filter[, "Cor.Diff"])
> jpeg("/home/ehsab/Gene_duplication/Test.jpg", width = 9, height=7, units = "in", res = 400)
> plot(TT.filter[, "Score"], TT.filter[, "Cor.Diff"])
> dev.off()
