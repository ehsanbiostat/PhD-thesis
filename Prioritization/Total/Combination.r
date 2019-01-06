
w=list()


load("1_text.RData")
load("2_text.RData")
load("3_text.RData")
load("4_text.RData")
load("5_text.RData")
load("6_text.RData")
load("7_text.RData")
for (i in 1:147){
 q=data.frame(alldeg, bet, clos, klein, apply(sim.jacard[,-i],1,sum), apply(sim.dic[,-i],1,sum), apply(sim.we[,-i],1,sum), apply(short[,-i],1,sum))
 colnames(q)[5:8]=c("sim.jacard","sim.dic","sim.we","short")
 w[i]=list(q)
}

save(w,file="text.RData")






load("1_AG.RData")
load("2_AG.RData")
load("3_AG.RData")
load("4_AG.RData")
load("5_AG.RData")
load("6_AG.RData")
load("7_AG.RData")
for (i in 1:147){
 q=data.frame(alldeg, bet, clos, klein, apply(sim.jacard[,-i],1,sum), apply(sim.dic[,-i],1,sum), apply(sim.we[,-i],1,sum), apply(short[,-i],1,sum))
 colnames(q)[5:8]=c("sim.jacard","sim.dic","sim.we","short")
 w[i]=list(q)
}

save(w,file="AG.RData")





load("1_Genie3.RData")
load("2_Genie3.RData")
load("3_Genie3.RData")
load("4_Genie3.RData")
load("5_Genie3.RData")
load("6_Genie3.RData")
load("7_Genie3.RData")
for (i in 1:147){
 q=data.frame(alldeg, bet, clos, klein, apply(sim.jacard[,-i],1,sum), apply(sim.dic[,-i],1,sum), apply(sim.we[,-i],1,sum), apply(short[,-i],1,sum))
 colnames(q)[5:8]=c("sim.jacard","sim.dic","sim.we","short")
 w[i]=list(q)
}

save(w,file="Genie3.RData")




load("1_Genemania.RData")
load("2_Genemania.RData")
load("3_Genemania.RData")
load("4_Genemania.RData")
load("5_Genemania.RData")
load("6_Genemania.RData")
load("7_Genemania.RData")
for (i in 1:147){
 q=data.frame(alldeg, bet, clos, klein, apply(sim.jacard[,-i],1,sum), apply(sim.dic[,-i],1,sum), apply(sim.we[,-i],1,sum), apply(short[,-i],1,sum))
 colnames(q)[5:8]=c("sim.jacard","sim.dic","sim.we","short")
 w[i]=list(q)
}

save(w,file="Genemania.RData")


load("1_Mutant.RData")
load("2_Mutant.RData")
load("3_Mutant.RData")
load("4_Mutant.RData")
load("5_Mutant.RData")
load("6_Mutant.RData")
load("7_Mutant.RData")
for (i in 1:147){
 q=data.frame(alldeg, bet, clos, klein, apply(sim.jacard[,-i],1,sum), apply(sim.dic[,-i],1,sum), apply(sim.we[,-i],1,sum), apply(short[,-i],1,sum))
 colnames(q)[5:8]=c("sim.jacard","sim.dic","sim.we","short")
 w[i]=list(q)
 print(i)
}
save(w,file="Mutant.RData")


load("1_PCC.RData")
load("2_PCC.RData")
load("3_PCC.RData")
load("4_PCC.RData")
load("5_PCC.RData")
load("6_PCC.RData")
load("7_PCC.RData")
for (i in 1:147){
 q=data.frame(alldeg, bet, clos, klein, apply(sim.jacard[,-i],1,sum), apply(sim.dic[,-i],1,sum), apply(sim.we[,-i],1,sum), apply(short[,-i],1,sum))
 colnames(q)[5:8]=c("sim.jacard","sim.dic","sim.we","short")
 w[i]=list(q)
}
save(w,file="PCC.RData")


load("1_PPI.RData")
load("2_PPI.RData")
load("3_PPI.RData")
load("4_PPI.RData")
load("5_PPI.RData")
load("6_PPI.RData")
load("7_PPI.RData")
for (i in 1:147){
 q=data.frame(alldeg, bet, clos, klein, apply(sim.jacard[,-i],1,sum), apply(sim.dic[,-i],1,sum), apply(sim.we[,-i],1,sum), apply(short[,-i],1,sum))
 colnames(q)[5:8]=c("sim.jacard","sim.dic","sim.we","short")
 w[i]=list(q)
}
save(w,file="PPI.RData")
