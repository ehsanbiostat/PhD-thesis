
# library(ConsensusClusterPlus)
source("/ngsprojects/iwt/ehsab/Gene_duplication/Scripts/ConsensusClusterPlus/ConsensusClusterPlus.r")
# --------------------------------------------------------------------------------------------------------------------------------------------------
# In this case, it occurs because the function that plots dendrogram objects (stats:::plotNode) is implemented using a recursive algorithm and the dendrogram object is deeply nested.
# Ultimately, the correct solution is to modify plotNode to use an iterative algorithm, which will prevent the recursion depth error from occuring.
# In the short term, it is possible to force stats:::plotNode to be run as interpreted code rather then byte-compiled code via a nasty hack.
# --------------------------------------------------------------------------------------------------------------------------------------------------

unByteCode <- function(fun)
    {
        FUN <- eval(parse(text=deparse(fun)))
        environment(FUN) <- environment(fun)
        FUN
    }

## Replace function definition inside of a locked environment **HACK** 
assignEdgewise <- function(name, env, value)
    {
        unlockBinding(name, env=env)
        assign( name, envir=env, value=value)
        lockBinding(name, env=env)
        invisible(value)
    }

## Replace byte-compiled function in a locked environment with an interpreted-code
## function
unByteCodeAssign <- function(fun)
    {
        name <- gsub('^.*::+','', deparse(substitute(fun)))
        FUN <- unByteCode(fun)
        retval <- assignEdgewise(name=name,
                                 env=environment(FUN),
                                 value=FUN
                                 )
        invisible(retval)
    }

## Use the above functions to convert stats:::plotNode to interpreted-code:
unByteCodeAssign(stats:::plotNode)

## Now raise the interpreted code recursion limit (you may need to adjust this,
##  decreasing if it uses to much memory, increasing if you get a recursion depth error ).
options(expressions=5e4)


# --------------------------------------------------------------------------------------------------------------------------------------------------
# 
# -------------------------------------------------------------------------------------------------------------------------------------------------- 



## Read condition names
# condition.names = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/condition_name.txt", header = T)
# NAME = names(table(condition.names))
# k = 1
# NAME.spe = NAME[k]

# ## Read text file with read.table function
# setwd(paste("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Merge_module/",NAME[k],"/Dist/", sep = ""))

# g = lapply(paste(1:19285,".txt", sep = ""), read.table)
# Dist = matrix(NA, length(unlist(g[1])), length(g))
# for (i in 1:length(g)){
# 	print(i)
# 	Dist[, i] =  as.matrix(unlist(g[i]))
# }

load("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Merge_module/abiotic/G1_dist.RData")
setwd("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Merge_module/abiotic/")
rcc = ConsensusClusterPlus(Dist, maxK = 100, reps = 20, pItem = 0.8, clusterAlg = "kmdist", plot = "png", writeTable = T, verbose = T)
resICL = calcICL(rcc,title="example", writeTable=T)