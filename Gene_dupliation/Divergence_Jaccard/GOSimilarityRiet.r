GOtable = select(org.At.tair.db, keys = keys(org.At.tair.db), columns = c("TAIR", "GOALL"))
GOtable.BP <- subset(GOtable, ONTOLOGYALL == "BP")
GOtable.MF <- subset(GOtable, ONTOLOGYALL == "MF")
GOtable.CC <- subset(GOtable, ONTOLOGYALL == "CC")
universe <- row.names(MA.data)
universe.BP <- intersect(universe, unique(GOtable.BP$TAIR)) # take the overlap of genes in the MA.data and measure
universe.MF <- intersect(universe, unique(GOtable.MF$TAIR)) # take the overlap of genes in the MA.data and measure
universe.CC <- intersect(universe, unique(GOtable.CC$TAIR)) # take the overlap of genes in the MA.data and measure
hgCutoff <- 0.0001
params <- new("GOHyperGParams",
geneIds=as.vector(module1G$mod1),
universeGeneIds=universe.BP,
annotation="org.At.tair.db",
ontology="BP",
pvalueCutoff=hgCutoff,
conditional=FALSE,
testDirection="over")
# Get the GO functions enriched in module 1
hgOver.BP.1 <- hyperGTest(params)
S.1 <- summary(hgOver.BP.1)
# Get the GO functions enriched in module 2
newParams <- params
geneIds(newParams) <- as.vector(module2G$mod2)
hgOver.BP.2 <- hyperGTest(newParams)
S.2 <- summary(hgOver.BP.2)
## Display the similarity of the GO-terms
S.1$module <- 1
S.2$module <- 2
GO_results <- rbind(S.1, S.2)
delRows <- c(which(GO_results$Size> 2000), which(GO_results$size < 30)) # Remove GO terms with very few or very many genes
GO_results <- GO_results[-delRows,]
GO_results$module[which(GO_results$GOBPID %in% intersect(S.1$GOBPID, S.2$GOBPID))]<- 3 # Modify the GO_results df such that module is 3 if both modules have the GO-terms
# Calculate the similarity of the GOterms
GO_sel <- subset(GOtable.BP, GOALL %in% GO_results$GOBPID, select = c("TAIR", "GOALL"))
GO_sel.tab <- table(GO_sel)
GO_sel.df <- as.data.frame.table(GO_sel.tab)
GO_gf.mat <- acast(GO_sel.df, GOALL ~ TAIR, value.var = "Freq")
#delCols <- which(colSums(GO_gf.mat) == 0)#remove columns with only zeros
#GO_gf.mat <- GO_gf.mat[,-delCols]
clus <- hclust(dist(GO_gf.mat, method = "binary"), method = "ward.D2")
moduleTips <- data.frame("module" = GO_results$module[match(clus$labels, GO_results$GOBPID)], "GO" = clus$labels)
moduleTips$module[which(moduleTips$GO %in% intersect(intersect(S.1$GOBPID, S.2$GOBPID), GO_results$GOBPID))]<- 3 #Indicate the GO-terms that refer to both clusters
colors <- c("#036564", "#EB6841", "#3366CC") # define the colors for the labels
colors <- c("#006400", "#00008B", "#BB0000") # define the colors for the labels
mapColors <- merge(moduleTips, data.frame("module" = c(1:3), "col" = colors), by = "module", all.x = TRUE)
colorsLabels <- mapColors$col[match(clus$labels, mapColors$GO)]
png(file.path(plotDir, paste("GODendro_", mod1, "_", mod2, ".png", sep = "")), units="px", width=600, height=1500, res = 300)
plot(as.phylo(clus), cex = 0.2, no.margin = TRUE, tip.color = as.vector(colorsLabels))
dev.off()
# TO DO: maybe add a bar chart to the graph which represents the -log(p) ==> significance of GO terms
write.table(GO_results[match(clus$labels[clus$order], GO_results$GOBPID),], file = file.path(plotDir, paste("GOResults_BP_", mod1, "_", mod2, ".txt", sep = "")), quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
