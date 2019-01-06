# Med.cond.score1 = apply(CondRank1, 2, median)
# Med.cond.score2 = apply(CondRank2, 2, median)
CondRank = rbind(c(1:dim(total_cond_scores)[1]), CondRank1, CondRank2)
colnames(CondRank) = colnames(data.exp.scale.med)
CondRank0 = CondRank[,order(CondRank[ORDER + 1, ], decreasing = T)]
CondRank = rbind(c(1:dim(total_cond_scores)[1]), CondRank2, CondRank1)
colnames(CondRank) = colnames(data.exp.scale.med)
CondRank01 = CondRank[,order(CondRank[ORDER + 1, ], decreasing = T)]
 
# Entire module heatmap
# jpeg(paste("/group/biocomp/users/ehsab/Gene_duplication/Results/SA/Local/Validation/Compendium_wise/Selected_V4/development/95/Plot/Module/Top10/",k,".jpeg", sep = ""), width = 11, height = 7, units = "in", res = 200, pointsize=8)
jpeg(paste("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Plot/Module/Top10/",k,".jpeg", sep = ""), width = 11, height = 7, units = "in", res = 200, pointsize=8)
   myImagePlot_sym(data.exp.scale[c(memb.comun, memb.comun.para),c(CondRank0[1,1:10],CondRank01[1,1:10])])
   dev.off()
 # Anchorpoints module heatmap
 # jpeg(paste("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Validation/Compendium_wise/development/Heatmap/Anchorpoints/",k,".jpeg", sep = ""), width = 11, height = 7, units = "in", res = 200, pointsize=8)
 # jpeg(paste("/group/biocomp/users/ehsab/Gene_duplication/Results/SA/Local/Validation/Compendium_wise/Selected_V4/development/95/Plot/Anchorpoint/Top10/",k,".jpeg", sep = ""), width = 11, height = 7, units = "in", res = 200, pointsize=8)
 jpeg(paste("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Plot/Anchorpoint/Top10/",k,".jpeg", sep = ""), width = 11, height = 7, units = "in", res = 200, pointsize=8)
 myImagePlot_sym(data.exp.scale[c(j, jpar),c(CondRank0[1,1:10],CondRank01[1,1:10])])
 dev.off()
 # ------------------------------------------------------------------
 # CORNET annotation
# Annot.j = foreach(i = 1:10, .combine = c) %do% {
# strsplit(substring(colnames(CondRank0)[1:10], 2), "[.]")[[i]][1]
# }
# Annot.jpar = foreach(i = 1:10, .combine = c) %do% {
# strsplit(substring(colnames(CondRank01)[1:10], 2), "[.]")[[i]][1]
# }
# ------------------------------------------------------------------
# data.frame(cbind(Annot.j, Annot.jpar), memb.comun[1], memb.comun.para[1])
}
top
Annot = foreach (N = 1:nrow(top), .combine = rbind) %dopar% {
   j = top[N, 1]
   jpar = top[N, 2]
  memb.comun = neighborhood(g,1,j, mode = c("out"))[[1]]
memb.comun.para = neighborhood(g,1,jpar, mode = c("out"))[[1]]
Inter1 = intersect(memb.comun.para, paralogous[paralogous$Ref %in% memb.comun, "Anch"])
Inter2 = intersect(memb.comun, paralogous[paralogous$Ref %in% memb.comun.para, "Anch"])
aa = paralogous[paralogous$Ref %in% Inter2, ]
bb = aa[aa$Anch %in% Inter1,]
colnames(bb)[1] = "References"
Inter_name = merge(ref, bb, by = "References")
colnames(Inter_name)[c(1,3)] = c("Ref", "References")
Inter_name = merge(ref, Inter_name, by = "References")
i = paste(ref[j, 1], ref[jpar, 1], top[N, 5], sep = "_")
# write.table(Inter_name, file = paste("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Local/Compendium_wise/",NAME.spe,"/Score99/Module_Anchorpoints_V2/Module",i,".txt", sep = ""), row.names = F, quote = F, sep = "\t")
# write.table(Inter_name, file = paste("/group/biocomp/users/ehsab/Gene_duplication/Results/SA/Local/Validation/Compendium_wise/Selected_V4/development/95/Genes/Anchorpoint/",i,".txt", sep = ""), row.names = F, quote = F, sep = "\t")
# write.table(ref[memb.comun, 1], file = paste("/group/biocomp/users/ehsab/Gene_duplication/Results/SA/Local/Validation/Compendium_wise/Selected_V4/development/95/Genes/Module/",ref[j, 1],".txt", sep = ""), col.names = F, row.names = F, quote = F, sep = "\t")
# write.table(ref[memb.comun.para, 1], file = paste("/group/biocomp/users/ehsab/Gene_duplication/Results/SA/Local/Validation/Compendium_wise/Selected_V4/development/95/Genes/Module/",ref[jpar, 1],".txt", sep = ""), col.names = F, row.names = F, quote = F, sep = "\t")
write.table(Inter_name, file = paste("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/Anchorpoint/",i,".txt", sep = ""), row.names = F, quote = F, sep = "\t")
write.table(ref[memb.comun, 1], file = paste("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/Module/",ref[j, 1],".txt", sep = ""), col.names = F, row.names = F, quote = F, sep = "\t")
write.table(ref[memb.comun.para, 1], file = paste("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/Module/",ref[jpar, 1],".txt", sep = ""), col.names = F, row.names = F, quote = F, sep = "\t")
colnames(Inter_name)[c(1,3)] = c("Ref", "References")
Result = merge(ref, Inter_name, by = "References")
# k = paste(ref[j, 1], ref[jpar, 1], Image_scoring[N, 5], sep = "_")
k = paste(ref[j, 1], ref[jpar, 1], top[N, 5], sep = "_")
j = Result$References #Gene.y
jpar = Result$Ref #Gene.x
CondRank1 = foreach(i = 1:length(memb.comun), .combine = rbind) %do%{
rank(total_cond_scores[, memb.comun[i]])
}
CondRank2 = foreach(i = 1:length(memb.comun.para), .combine = rbind) %do%{
rank(total_cond_scores[, memb.comun.para[i]])
}
ORDER = 1 
# Med.cond.score1 = apply(CondRank1, 2, median)
# Med.cond.score2 = apply(CondRank2, 2, median)
CondRank = rbind(c(1:dim(total_cond_scores)[1]), CondRank1, CondRank2)
colnames(CondRank) = colnames(data.exp.scale.med)
CondRank0 = CondRank[,order(CondRank[ORDER + 1, ], decreasing = T)]
CondRank = rbind(c(1:dim(total_cond_scores)[1]), CondRank2, CondRank1)
colnames(CondRank) = colnames(data.exp.scale.med)
CondRank01 = CondRank[,order(CondRank[ORDER + 1, ], decreasing = T)]
 
# Entire module heatmap
# jpeg(paste("/group/biocomp/users/ehsab/Gene_duplication/Results/SA/Local/Validation/Compendium_wise/Selected_V4/development/95/Plot/Module/Top10/",k,".jpeg", sep = ""), width = 11, height = 7, units = "in", res = 200, pointsize=8)
jpeg(paste("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Plot/Module/Top10/",k,".jpeg", sep = ""), width = 11, height = 7, units = "in", res = 200, pointsize=8)
   myImagePlot_sym(data.exp.scale[c(memb.comun, memb.comun.para),c(CondRank0[1,1:10],CondRank01[1,1:10])])
   dev.off()
 # Anchorpoints module heatmap
 # jpeg(paste("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Validation/Compendium_wise/development/Heatmap/Anchorpoints/",k,".jpeg", sep = ""), width = 11, height = 7, units = "in", res = 200, pointsize=8)
 # jpeg(paste("/group/biocomp/users/ehsab/Gene_duplication/Results/SA/Local/Validation/Compendium_wise/Selected_V4/development/95/Plot/Anchorpoint/Top10/",k,".jpeg", sep = ""), width = 11, height = 7, units = "in", res = 200, pointsize=8)
 jpeg(paste("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Plot/Anchorpoint/Top10/",k,".jpeg", sep = ""), width = 11, height = 7, units = "in", res = 200, pointsize=8)
 myImagePlot_sym(data.exp.scale[c(j, jpar),c(CondRank0[1,1:10],CondRank01[1,1:10])])
 dev.off()
 # ------------------------------------------------------------------
 # CORNET annotation
# Annot.j = foreach(i = 1:10, .combine = c) %do% {
# strsplit(substring(colnames(CondRank0)[1:10], 2), "[.]")[[i]][1]
# }
# Annot.jpar = foreach(i = 1:10, .combine = c) %do% {
# strsplit(substring(colnames(CondRank01)[1:10], 2), "[.]")[[i]][1]
# }
# ------------------------------------------------------------------
# data.frame(cbind(Annot.j, Annot.jpar), memb.comun[1], memb.comun.para[1])
}
top
Annot = foreach (N = 1:nrow(top), .combine = rbind) %dopar% {
   j = top[N, 1]
   jpar = top[N, 2]
  memb.comun = neighborhood(g,1,j, mode = c("out"))[[1]]
memb.comun.para = neighborhood(g,1,jpar, mode = c("out"))[[1]]
Inter1 = intersect(memb.comun.para, paralogous[paralogous$Ref %in% memb.comun, "Anch"])
Inter2 = intersect(memb.comun, paralogous[paralogous$Ref %in% memb.comun.para, "Anch"])
aa = paralogous[paralogous$Ref %in% Inter2, ]
bb = aa[aa$Anch %in% Inter1,]
colnames(bb)[1] = "References"
Inter_name = merge(ref, bb, by = "References")
colnames(Inter_name)[c(1,3)] = c("Ref", "References")
Inter_name = merge(ref, Inter_name, by = "References")
i = paste(ref[j, 1], ref[jpar, 1], top[N, 5], top[N, 6], sep = "_")
# write.table(Inter_name, file = paste("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Local/Compendium_wise/",NAME.spe,"/Score99/Module_Anchorpoints_V2/Module",i,".txt", sep = ""), row.names = F, quote = F, sep = "\t")
# write.table(Inter_name, file = paste("/group/biocomp/users/ehsab/Gene_duplication/Results/SA/Local/Validation/Compendium_wise/Selected_V4/development/95/Genes/Anchorpoint/",i,".txt", sep = ""), row.names = F, quote = F, sep = "\t")
# write.table(ref[memb.comun, 1], file = paste("/group/biocomp/users/ehsab/Gene_duplication/Results/SA/Local/Validation/Compendium_wise/Selected_V4/development/95/Genes/Module/",ref[j, 1],".txt", sep = ""), col.names = F, row.names = F, quote = F, sep = "\t")
# write.table(ref[memb.comun.para, 1], file = paste("/group/biocomp/users/ehsab/Gene_duplication/Results/SA/Local/Validation/Compendium_wise/Selected_V4/development/95/Genes/Module/",ref[jpar, 1],".txt", sep = ""), col.names = F, row.names = F, quote = F, sep = "\t")
write.table(Inter_name, file = paste("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/Anchorpoint/",i,".txt", sep = ""), row.names = F, quote = F, sep = "\t")
write.table(ref[memb.comun, 1], file = paste("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/Module/",ref[j, 1],"_",top[N, 6],".txt", sep = ""), col.names = F, row.names = F, quote = F, sep = "\t")
write.table(ref[memb.comun.para, 1], file = paste("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/Module/",ref[jpar, 1],"_",top[N, 6]".txt", sep = ""), col.names = F, row.names = F, quote = F, sep = "\t")
colnames(Inter_name)[c(1,3)] = c("Ref", "References")
Result = merge(ref, Inter_name, by = "References")
# k = paste(ref[j, 1], ref[jpar, 1], Image_scoring[N, 5], sep = "_")
k = paste(ref[j, 1], ref[jpar, 1], top[N, 5], top[N, 6], sep = "_")
j = Result$References #Gene.y
jpar = Result$Ref #Gene.x
CondRank1 = foreach(i = 1:length(memb.comun), .combine = rbind) %do%{
rank(total_cond_scores[, memb.comun[i]])
}
CondRank2 = foreach(i = 1:length(memb.comun.para), .combine = rbind) %do%{
rank(total_cond_scores[, memb.comun.para[i]])
}
ORDER = 1 
# Med.cond.score1 = apply(CondRank1, 2, median)
# Med.cond.score2 = apply(CondRank2, 2, median)
CondRank = rbind(c(1:dim(total_cond_scores)[1]), CondRank1, CondRank2)
colnames(CondRank) = colnames(data.exp.scale.med)
CondRank0 = CondRank[,order(CondRank[ORDER + 1, ], decreasing = T)]
Annot = foreach (N = 1:nrow(top), .combine = rbind) %dopar% {
   j = top[N, 1]
   jpar = top[N, 2]
  memb.comun = neighborhood(g,1,j, mode = c("out"))[[1]]
memb.comun.para = neighborhood(g,1,jpar, mode = c("out"))[[1]]
Inter1 = intersect(memb.comun.para, paralogous[paralogous$Ref %in% memb.comun, "Anch"])
Inter2 = intersect(memb.comun, paralogous[paralogous$Ref %in% memb.comun.para, "Anch"])
aa = paralogous[paralogous$Ref %in% Inter2, ]
bb = aa[aa$Anch %in% Inter1,]
colnames(bb)[1] = "References"
Inter_name = merge(ref, bb, by = "References")
colnames(Inter_name)[c(1,3)] = c("Ref", "References")
Inter_name = merge(ref, Inter_name, by = "References")
i = paste(ref[j, 1], ref[jpar, 1], top[N, 5], top[N, 6], sep = "_")
# write.table(Inter_name, file = paste("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Local/Compendium_wise/",NAME.spe,"/Score99/Module_Anchorpoints_V2/Module",i,".txt", sep = ""), row.names = F, quote = F, sep = "\t")
# write.table(Inter_name, file = paste("/group/biocomp/users/ehsab/Gene_duplication/Results/SA/Local/Validation/Compendium_wise/Selected_V4/development/95/Genes/Anchorpoint/",i,".txt", sep = ""), row.names = F, quote = F, sep = "\t")
# write.table(ref[memb.comun, 1], file = paste("/group/biocomp/users/ehsab/Gene_duplication/Results/SA/Local/Validation/Compendium_wise/Selected_V4/development/95/Genes/Module/",ref[j, 1],".txt", sep = ""), col.names = F, row.names = F, quote = F, sep = "\t")
# write.table(ref[memb.comun.para, 1], file = paste("/group/biocomp/users/ehsab/Gene_duplication/Results/SA/Local/Validation/Compendium_wise/Selected_V4/development/95/Genes/Module/",ref[jpar, 1],".txt", sep = ""), col.names = F, row.names = F, quote = F, sep = "\t")
write.table(Inter_name, file = paste("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/Anchorpoint/",i,".txt", sep = ""), row.names = F, quote = F, sep = "\t")
write.table(ref[memb.comun, 1], file = paste("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/Module/",ref[j, 1],"_",top[N, 6],".txt", sep = ""), col.names = F, row.names = F, quote = F, sep = "\t")
write.table(ref[memb.comun.para, 1], file = paste("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/Module/",ref[jpar, 1],"_",top[N, 6],".txt", sep = ""), col.names = F, row.names = F, quote = F, sep = "\t")
colnames(Inter_name)[c(1,3)] = c("Ref", "References")
Result = merge(ref, Inter_name, by = "References")
# k = paste(ref[j, 1], ref[jpar, 1], Image_scoring[N, 5], sep = "_")
k = paste(ref[j, 1], ref[jpar, 1], top[N, 5], top[N, 6], sep = "_")
j = Result$References #Gene.y
jpar = Result$Ref #Gene.x
CondRank1 = foreach(i = 1:length(memb.comun), .combine = rbind) %do%{
rank(total_cond_scores[, memb.comun[i]])
}
CondRank2 = foreach(i = 1:length(memb.comun.para), .combine = rbind) %do%{
rank(total_cond_scores[, memb.comun.para[i]])
}
ORDER = 1 
# Med.cond.score1 = apply(CondRank1, 2, median)
# Med.cond.score2 = apply(CondRank2, 2, median)
CondRank = rbind(c(1:dim(total_cond_scores)[1]), CondRank1, CondRank2)
colnames(CondRank) = colnames(data.exp.scale.med)
CondRank0 = CondRank[,order(CondRank[ORDER + 1, ], decreasing = T)]
CondRank = rbind(c(1:dim(total_cond_scores)[1]), CondRank2, CondRank1)
colnames(CondRank) = colnames(data.exp.scale.med)
CondRank01 = CondRank[,order(CondRank[ORDER + 1, ], decreasing = T)]
 
# Entire module heatmap
# jpeg(paste("/group/biocomp/users/ehsab/Gene_duplication/Results/SA/Local/Validation/Compendium_wise/Selected_V4/development/95/Plot/Module/Top10/",k,".jpeg", sep = ""), width = 11, height = 7, units = "in", res = 200, pointsize=8)
jpeg(paste("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Plot/Module/Top10/",k,".jpeg", sep = ""), width = 11, height = 7, units = "in", res = 200, pointsize=8)
   myImagePlot_sym(data.exp.scale[c(memb.comun, memb.comun.para),c(CondRank0[1,1:10],CondRank01[1,1:10])])
   dev.off()
 # Anchorpoints module heatmap
 # jpeg(paste("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Validation/Compendium_wise/development/Heatmap/Anchorpoints/",k,".jpeg", sep = ""), width = 11, height = 7, units = "in", res = 200, pointsize=8)
 # jpeg(paste("/group/biocomp/users/ehsab/Gene_duplication/Results/SA/Local/Validation/Compendium_wise/Selected_V4/development/95/Plot/Anchorpoint/Top10/",k,".jpeg", sep = ""), width = 11, height = 7, units = "in", res = 200, pointsize=8)
 jpeg(paste("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Plot/Anchorpoint/Top10/",k,".jpeg", sep = ""), width = 11, height = 7, units = "in", res = 200, pointsize=8)
 myImagePlot_sym(data.exp.scale[c(j, jpar),c(CondRank0[1,1:10],CondRank01[1,1:10])])
 dev.off()
 # ------------------------------------------------------------------
 # CORNET annotation
# Annot.j = foreach(i = 1:10, .combine = c) %do% {
# strsplit(substring(colnames(CondRank0)[1:10], 2), "[.]")[[i]][1]
# }
# Annot.jpar = foreach(i = 1:10, .combine = c) %do% {
# strsplit(substring(colnames(CondRank01)[1:10], 2), "[.]")[[i]][1]
# }
# ------------------------------------------------------------------
# data.frame(cbind(Annot.j, Annot.jpar), memb.comun[1], memb.comun.para[1])
}
top = a[a$Group == 90,]
top
load(paste("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Gene_score/Compendium_wise/",NAME[k],"/Graph_V2_90.RData", sep = ""))
k = 6
load(paste("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Gene_score/Compendium_wise/",NAME[k],"/Graph_V2_90.RData", sep = ""))
g
Annot = foreach (N = 1:nrow(top), .combine = rbind) %dopar% {
   j = top[N, 1]
   jpar = top[N, 2]
  memb.comun = neighborhood(g,1,j, mode = c("out"))[[1]]
memb.comun.para = neighborhood(g,1,jpar, mode = c("out"))[[1]]
Inter1 = intersect(memb.comun.para, paralogous[paralogous$Ref %in% memb.comun, "Anch"])
Inter2 = intersect(memb.comun, paralogous[paralogous$Ref %in% memb.comun.para, "Anch"])
aa = paralogous[paralogous$Ref %in% Inter2, ]
bb = aa[aa$Anch %in% Inter1,]
colnames(bb)[1] = "References"
Inter_name = merge(ref, bb, by = "References")
colnames(Inter_name)[c(1,3)] = c("Ref", "References")
Inter_name = merge(ref, Inter_name, by = "References")
i = paste(ref[j, 1], ref[jpar, 1], top[N, 5], top[N, 6], sep = "_")
# write.table(Inter_name, file = paste("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Local/Compendium_wise/",NAME.spe,"/Score99/Module_Anchorpoints_V2/Module",i,".txt", sep = ""), row.names = F, quote = F, sep = "\t")
# write.table(Inter_name, file = paste("/group/biocomp/users/ehsab/Gene_duplication/Results/SA/Local/Validation/Compendium_wise/Selected_V4/development/95/Genes/Anchorpoint/",i,".txt", sep = ""), row.names = F, quote = F, sep = "\t")
# write.table(ref[memb.comun, 1], file = paste("/group/biocomp/users/ehsab/Gene_duplication/Results/SA/Local/Validation/Compendium_wise/Selected_V4/development/95/Genes/Module/",ref[j, 1],".txt", sep = ""), col.names = F, row.names = F, quote = F, sep = "\t")
# write.table(ref[memb.comun.para, 1], file = paste("/group/biocomp/users/ehsab/Gene_duplication/Results/SA/Local/Validation/Compendium_wise/Selected_V4/development/95/Genes/Module/",ref[jpar, 1],".txt", sep = ""), col.names = F, row.names = F, quote = F, sep = "\t")
write.table(Inter_name, file = paste("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/Anchorpoint/",i,".txt", sep = ""), row.names = F, quote = F, sep = "\t")
write.table(ref[memb.comun, 1], file = paste("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/Module/",ref[j, 1],"_",top[N, 6],".txt", sep = ""), col.names = F, row.names = F, quote = F, sep = "\t")
write.table(ref[memb.comun.para, 1], file = paste("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/Module/",ref[jpar, 1],"_",top[N, 6],".txt", sep = ""), col.names = F, row.names = F, quote = F, sep = "\t")
colnames(Inter_name)[c(1,3)] = c("Ref", "References")
Result = merge(ref, Inter_name, by = "References")
# k = paste(ref[j, 1], ref[jpar, 1], Image_scoring[N, 5], sep = "_")
k = paste(ref[j, 1], ref[jpar, 1], top[N, 5], top[N, 6], sep = "_")
j = Result$References #Gene.y
jpar = Result$Ref #Gene.x
CondRank1 = foreach(i = 1:length(memb.comun), .combine = rbind) %do%{
rank(total_cond_scores[, memb.comun[i]])
}
CondRank2 = foreach(i = 1:length(memb.comun.para), .combine = rbind) %do%{
rank(total_cond_scores[, memb.comun.para[i]])
}
ORDER = 1 
# Med.cond.score1 = apply(CondRank1, 2, median)
# Med.cond.score2 = apply(CondRank2, 2, median)
CondRank = rbind(c(1:dim(total_cond_scores)[1]), CondRank1, CondRank2)
colnames(CondRank) = colnames(data.exp.scale.med)
CondRank0 = CondRank[,order(CondRank[ORDER + 1, ], decreasing = T)]
CondRank = rbind(c(1:dim(total_cond_scores)[1]), CondRank2, CondRank1)
colnames(CondRank) = colnames(data.exp.scale.med)
CondRank01 = CondRank[,order(CondRank[ORDER + 1, ], decreasing = T)]
 
# Entire module heatmap
# jpeg(paste("/group/biocomp/users/ehsab/Gene_duplication/Results/SA/Local/Validation/Compendium_wise/Selected_V4/development/95/Plot/Module/Top10/",k,".jpeg", sep = ""), width = 11, height = 7, units = "in", res = 200, pointsize=8)
jpeg(paste("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Plot/Module/Top10/",k,".jpeg", sep = ""), width = 11, height = 7, units = "in", res = 200, pointsize=8)
   myImagePlot_sym(data.exp.scale[c(memb.comun, memb.comun.para),c(CondRank0[1,1:10],CondRank01[1,1:10])])
   dev.off()
 # Anchorpoints module heatmap
 # jpeg(paste("/scratch/iw
 # jpeg(paste("/group/biocomp/users/ehsab/Gene_duplicati
 
 
 dev.off()
 # ------------------------------------------------------------------
 # CORNET annotation
# Annot.j = foreach(i = 1:10, .combine = c) %do% {
# strsplit(substring(colnames(CondRank0)[1:10], 2), "[.]")[[
# }
# Annot.jpar = foreach(i = 1:10, .combine = c) %do% {
# strsplit(substring(colnames(CondRank01)[1:10]
# }
# ------------------------------------------------------------------
# data.frame(cbind(Annot.j, Annot.j
}
}
}
Annot = foreach (N = 1:nrow(top), .combine = rbind) %dopar% {
   j = top[N, 1]
   jpar = top[N, 2]
  memb.comun = neighborhood(g,1,j, mode = c("out"))[[1]]
memb.comun.para = neighborhood(g,1,jpar, mode = c("out"))[[1]]
Inter1 = intersect(memb.comun.para, paralogous[paralogous$Ref %in% memb.comun, "Anch"])
Inter2 = intersect(memb.comun, paralogous[paralogous$Ref %in% memb.comun.para, "Anch"])
aa = paralogous[paralogous$Ref %in% Inter2, ]
bb = aa[aa$Anch %in% Inter1,]
colnames(bb)[1] = "References"
Inter_name = merge(ref, bb, by = "References")
colnames(Inter_name)[c(1,3)] = c("Ref", "References")
Inter_name = merge(ref, Inter_name, by = "References")
i = paste(ref[j, 1], ref[jpar, 1], top[N, 5], top[N, 6], sep = "_")
# write.table(Inter_name, file = paste("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Local/Compendium_wise/",NAME.spe,"/Score99/Module_Anchorpoints_V2/Module",i,".txt", sep = ""), row.names = F, quote = F, sep = "\t")
# write.table(Inter_name, file = paste("/group/biocomp/users/ehsab/Gene_duplication/Results/SA/Local/Validation/Compendium_wise/Selected_V4/development/95/Genes/Anchorpoint/",i,".txt", sep = ""), row.names = F, quote = F, sep = "\t")
# write.table(ref[memb.comun, 1], file = paste("/group/biocomp/users/ehsab/Gene_duplication/Results/SA/Local/Validation/Compendium_wise/Selected_V4/development/95/Genes/Module/",ref[j, 1],".txt", sep = ""), col.names = F, row.names = F, quote = F, sep = "\t")
# write.table(ref[memb.comun.para, 1], file = paste("/group/biocomp/users/ehsab/Gene_duplication/Results/SA/Local/Validation/Compendium_wise/Selected_V4/development/95/Genes/Module/",ref[jpar, 1],".txt", sep = ""), col.names = F, row.names = F, quote = F, sep = "\t")
write.table(Inter_name, file = paste("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/Anchorpoint/",i,".txt", sep = ""), row.names = F, quote = F, sep = "\t")
write.table(ref[memb.comun, 1], file = paste("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/Module/",ref[j, 1],"_",top[N, 6],".txt", sep = ""), col.names = F, row.names = F, quote = F, sep = "\t")
write.table(ref[memb.comun.para, 1], file = paste("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/Module/",ref[jpar, 1],"_",top[N, 6],".txt", sep = ""), col.names = F, row.names = F, quote = F, sep = "\t")
colnames(Inter_name)[c(1,3)] = c("Ref", "References")
Result = merge(ref, Inter_name, by = "References")
# k = paste(ref[j, 1], ref[jpar, 1], Image_scoring[N, 5], sep = "_")
k = paste(ref[j, 1], ref[jpar, 1], top[N, 5], top[N, 6], sep = "_")
j = Result$References #Gene.y
jpar = Result$Ref #Gene.x
CondRank1 = foreach(i = 1:length(memb.comun), .combine = rbind) %do%{
rank(total_cond_scores[, memb.comun[i]])
}
CondRank2 = foreach(i = 1:length(memb.comun.para), .combine = rbind) %do%{
rank(total_cond_scores[, memb.comun.para[i]])
}
ORDER = 1 
# Med.cond.score1 = apply(CondRank1, 2, median)
# Med.cond.score2 = apply(CondRank2, 2, median)
CondRank = rbind(c(1:dim(total_cond_scores)[1]), CondRank1, CondRank2)
colnames(CondRank) = colnames(data.exp.scale.med)
CondRank0 = CondRank[,order(CondRank[ORDER + 1, ], decreasing = T)]
CondRank = rbind(c(1:dim(total_cond_scores)[1]), CondRank2, CondRank1)
colnames(CondRank) = colnames(data.exp.scale.med)
CondRank01 = CondRank[,order(CondRank[ORDER + 1, ], decreasing = T)]
 
# Entire module heatmap
# jpeg(paste("/group/biocomp/users/ehsab/Gene_duplication/Results/SA/Local/Validation/Compendium_wise/Selected_V4/development/95/Plot/Module/Top10/",k,".jpeg", sep = ""), width = 11, height = 7, units = "in", res = 200, pointsize=8)
jpeg(paste("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Plot/Module/Top10/",k,".jpeg", sep = ""), width = 11, height = 7, units = "in", res = 200, pointsize=8)
   myImagePlot_sym(data.exp.scale[c(memb.comun, memb.comun.para),c(CondRank0[1,1:10],CondRank01[1,1:10])])
   dev.off()
 # Anchorpoints module heatmap
 # jpeg(paste("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Validation/Compendium_wise/development/Heatmap/Anchorpoints/",k,".jpeg", sep = ""), width = 11, height = 7, units = "in", res = 200, pointsize=8)
 # jpeg(paste("/group/biocomp/users/ehsab/Gene_duplication/Results/SA/Local/Validation/Compendium_wise/Selected_V4/development/95/Plot/Anchorpoint/Top10/",k,".jpeg", sep = ""), width = 11, height = 7, units = "in", res = 200, pointsize=8)
 jpeg(paste("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Plot/Anchorpoint/Top10/",k,".jpeg", sep = ""), width = 11, height = 7, units = "in", res = 200, pointsize=8)
 myImagePlot_sym(data.exp.scale[c(j, jpar),c(CondRank0[1,1:10],CondRank01[1,1:10])])
 dev.off()
 # ------------------------------------------------------------------
 # CORNET annotation
# Annot.j = foreach(i = 1:10, .combine = c) %do% {
# strsplit(substring(colnames(CondRank0)[1:10], 2), "[.]")[[i]][1]
# }
# Annot.jpar = foreach(i = 1:10, .combine = c) %do% {
# strsplit(substring(colnames(CondRank01)[1:10], 2), "[.]")[[i]][1]
# }
# ------------------------------------------------------------------
# data.frame(cbind(Annot.j, Annot.jpar), memb.comun[1], memb.comun.para[1])
}
top
N
j = top[N, 1]
   jpar = top[N, 2]
  memb.comun = neighborhood(g,1,j, mode = c("out"))[[1]]
memb.comun.para = neighborhood(g,1,jpar, mode = c("out"))[[1]]
N = 3
Gene
Gene = unique(c(top$Gene, top$Paralogous))
Gene
memb.comun = neighborhood(g,1,Gene, mode = c("out"))[[1]]
memb.comun
memb.comun = neighborhood(g,1,Gene, mode = c("out"))
memb.comun
Jacc = matrix(NA, length(Gene), length(Gene))
Jacc
colnames(Jacc) = rownames(Jacc) = Gene
for (i in 1:length(Gene)) { 
memb.comun = neighborhood(g,1,Gene, mode = c("out"))[[i]]
for (j in 1:length(Gene)) { 
memb.comun.para = neighborhood(g,1,Gene, mode = c("out"))[[j]]
Jacc[i, j] = length(intersect(memb.comun, memb.comun.para))/length(union(memb.comun, memb.comun.para))
}
}
plot(hclust(as.dist(1-Jacc)))
library(corrplot)
corrplot(Jacc)
str(corrplot)
corrplot(Jacc, method = c("color"), is.corr = F, order = c("AOE"))
corrplot(Jacc, method = c("number"), is.corr = F, order = c("AOE"))
Gene
memb.comun = unlist(neighborhood(g,1,Gene[c(2,3,4)], mode = c("out")))
memb.comun
sort(table(memb.comun))
Jacc
plot(hclust(as.dist(1-Jacc)))
library(venneuler)
install.packages("venneuler")
library(venneuler)
venneuler(Jacc)
plot(venneuler(Jacc))
Jacc
memb.comun = unlist(neighborhood(g,1,Gene, mode = c("out")))
GeneSet1 = Gene[c(1,6,5)]
universe = unlist(neighborhood(g,1,GeneSet1, mode = c("out")))
Counts <- matrix(NA, nrow=length(universe), ncol=length(GeneSet1))
dim(Counts)
for (i in 1:length(universe)) {
print(i)
for (j in 1:length(GeneSet1)) {
Counts[i,j] <- universe[i] %in% neighborhood(g,1,GeneSet1, mode = c("out"))[[j]]
}
}
GeneSet1.name = Gene.name[c(1,6,5)]
colnames(Counts) <- GeneSet1.name
library(limma)
cols<-c("Red", "Green", "Blue")
vennDiagram(vennCounts(Counts), circle.col=cols)
vennCounts(Counts[, c(1, 5, 6)])
class(vennCounts(Counts[, c(1, 5, 6)]))
names(vennCounts(Counts[, c(1, 5, 6)]))
Gene.name = ref[Gene, 2]
Gene.name
Gene.name = ref[Gene, 1]
Gene.name
colnames(Jacc) = rownames(Jacc) = Gene.name
plot(hclust(as.dist(1-Jacc)))
jpeg(paste("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/OverlappClustering95.jpeg", sep = ""), width = 11, height = 7, units = "in", res = 100, pointsize=8)
plot(hclust(as.dist(1-Jacc)))
dev.off()
jpeg(paste("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/OverlappClustering95.jpeg", sep = ""), width = 11, height = 7, units = "in", res = 200, pointsize=8)
plot(hclust(as.dist(1-Jacc)))
dev.off()
Gene.name
ref[12450,]
ref[2858,]
colnames(Counts) <- Gene.name
jpeg(paste("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/OverlappVenn95_1.jpeg", sep = ""), width = 11, height = 7, units = "in", res = 200, pointsize=8)
vennDiagram(vennCounts(Counts[, c(1, 5, 6)]), circle.col=cols)
dev.off()
Gene.name
vennDiagram(vennCounts(Counts[, c(2, 3, 4)]), circle.col=cols)
jpeg(paste("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/OverlappVenn95_2.jpeg", sep = ""), width = 11, height = 7, units = "in", res = 200, pointsize=8)
vennDiagram(vennCounts(Counts[, c(2, 3, 4)]), circle.col=cols)
dev.off()
vennDiagram(vennCounts(Counts[, c(7, 8)]), circle.col=cols)
top
g
Gene.name
Gene
degree(g, 9603)
degree(g, 9603, mode = c("out"))
degree(g, 2858, mode = c("out"))
vennCounts(Counts[, c(7, 8)])
vennCounts(Counts)
vennCounts(Counts[, c(7, 8)])[,1:2]
vennCounts(Counts[, c(7, 8)])[,1:5]
vennCounts(Counts)[,1:2]
vennCounts(Counts)[,7:8]
vennCounts(Counts)[,7:9]
head(vennCounts(Counts))
vennCounts(Counts)
Jacc.count = matrix(NA, length(Gene), length(Gene))
for (i in 1:length(Gene)) { 
memb.comun = neighborhood(g,1,Gene, mode = c("out"))[[i]]
for (j in 1:length(Gene)) { 
Jacc.count[i, j] = length(intersect(memb.comun, memb.comun.para))
}
}
for (i in 1:length(Gene)) { 
memb.comun = neighborhood(g,1,Gene, mode = c("out"))[[i]]
for (j in 1:length(Gene)) { 
memb.comun.para = neighborhood(g,1,Gene, mode = c("out"))[[j]]
Jacc.count[i, j] = length(intersect(memb.comun, memb.comun.para))
}
}
Jacc.count
colnames(Jacc.count) = rownames(Jacc.count) = Gene.name
Jacc.count
1300 + 1241 + 50 +393
savehistory("/ngsprojects//iwt/ehsab/Gene_duplication/Scripts/SA/History/ComparingDifferentThresholdPermutationSA.r")
