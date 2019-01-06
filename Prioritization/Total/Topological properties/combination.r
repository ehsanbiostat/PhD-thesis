load("~/prioritization/Leave-one-out/New/Total_topology/Bet.RData")
load("~/prioritization/Leave-one-out/New/Total_topology/closeness.RData")
load("~/prioritization/Leave-one-out/New/Total_topology/degree_klein.RData")
load("~/prioritization/Leave-one-out/New/Total_topology/dic.RData")
load("~/prioritization/Leave-one-out/New/Total_topology/jacard.RData")
load("~/prioritization/Leave-one-out/New/Total_topology/short.RData")
load("~/prioritization/Leave-one-out/New/Total_topology/we.RData")
a=data.frame(alldeg,decregiyg,cregiyg_sh,bet,clos,klein)

# combined netwrok
b=data.frame(a,apply(sim.jacard,1,sum),apply(sim.dic,1,sum),apply(sim.we,1,sum),apply(short,1,sum))
colnames(b)[7:10]=c("jaccard","sim.dic","sim.we","short")
write.table(b,file="~/prioritization/Leave-one-out/New/Total_topology/global.txt",sep="\t",row.names=F,col.names=T,quote=F)

# Leave-one-out datasets
d=list()
for (i in 1:147){
print(i)
a1=data.frame()
a1=cbind(a,apply(sim.jacard[,-i],1,sum))
a1=cbind(a1,apply(sim.dic[,-i],1,sum))
a1=cbind(a1,apply(sim.we[,-i],1,sum))
a1=cbind(a1,apply(short[,-i],1,sum))
colnames(a1)[7:10]=c("jaccard","sim.dic","sim.we","short")
d[i]=list(a1)
}
save(d, file="~/prioritization/Leave-one-out/New/Total_topology/Leave-one-out.RData")