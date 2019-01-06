
library(igraph)
library("isa2")
library("foreach")
library("doMC")
options(width = 200)

Score = read.delim("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Compendium_wise/development/PollenRoot90.txt")
Score95 = read.delim("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Compendium_wise/development/PollenRoot95.txt")
Score99 = read.delim("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Compendium_wise/development/PollenRoot99.txt")

PollenLeaf = intersect(Score95[Score95$PollenLeafPvalue < 0.01 & Score95$RootLeaf == 0, "Grape"], Score[Score$PollenLeafPvalue < 0.01 & Score$RootLeaf == 0, "Grape"])

Score[Score$Grape %in% PollenLeaf,]



intersect(Score[Score$RootLeafPvalue < 0.01, "Grape"], intersect(Score95[Score95$RootLeafPvalue < 0.01, "Grape"], Score99[Score99$RootLeafPvalue < 0.01, "Grape"]))
intersect(Score[Score$RootLeafPvalue < 0.001, "Grape"], intersect(Score95[Score95$RootLeafPvalue < 0.01, "Grape"], Score99[Score99$RootLeafPvalue < 0.01, "Grape"]))
intersect(Score[Score$RootLeafPvalue < 0.001, "Grape"], intersect(Score95[Score95$RootLeafPvalue < 0.001, "Grape"], Score99[Score99$RootLeafPvalue < 0.01, "Grape"]))
intersect(Score[Score$RootLeafPvalue < 0.001, "Grape"], intersect(Score95[Score95$RootLeafPvalue < 0.01, "Grape"], Score99[Score99$RootLeafPvalue < 0.01, "Grape"]))
Score[Score$Grape %in% c(7037, 7907, 14972),]
library(igraph)
load("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Graph_V2_90.RData")
intersect(neighborhood(g, 1, 7037, mode = c("out"))[[1]], intersect(neighborhood(g, 1, 7907, mode = c("out"))[[1]], neighborhood(g, 1, 14972, mode = c("out"))[[1]]))
options(width = 200)
intersect(neighborhood(g, 1, 7037, mode = c("out"))[[1]], intersect(neighborhood(g, 1, 7907, mode = c("out"))[[1]], neighborhood(g, 1, 14972, mode = c("out"))[[1]]))
neighborhood(g, 1, 7907, mode = c("out"))[[1]]
intersect(neighborhood(g, 1, 7037, mode = c("out"))[[1]], intersect(neighborhood(g, 1, 7907, mode = c("out"))[[1]], neighborhood(g, 1, 14972, mode = c("out"))[[1]]))
Score[Score$Grape %in% c(7037, 7907, 14972),]
PollenRootLeafclean = read.delim("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Orthology/Grape-Arabidopsis/PollenRootLeafclean.txt")
Triangle = read.delim("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Orthology/Grape-Arabidopsis/Triangle.txt")
Triangle
PollenRootLeafclean
Score[Score$Grape %in% c(7037, 7907, 14972),]
Score[Score$Grape %in% c(13166, 14526, 16881),]
Score[Score$Grape %in% c(7037, 7907, 14972),]
Score[Score$RootLeafPvalue < 0.0001, ]
Score[Score$PollenLeafPvalue < 0.0001, ]
savehistory("/home/ehsab/Gene_duplication/Scripts/history1.r")
