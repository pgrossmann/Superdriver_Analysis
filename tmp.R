load("/links/grid/scratch/grossman/rdata/B_mutator1000driver0.01factor2ratenormal0.00000001.rda")
source("analysis/scripts/analysis_functions.R")
if (!file.exists("tmp.rdata")) {
tmp <- getTime(B_mutator1000driver0.01factor2ratenormal0.00000001,c(3,2),whichValue=c("number.drivers","number.superdrivers"))
save(tmp,file="tmp.rdata")
} else load("tmp.rdata")
if (!file.exists("B_switched.rda")) {
B_switched <- switchPops2(B_mutator1000driver0.01factor2ratenormal0.00000001)
save(B_switched,file="B_switched.rda")
} else load("B_switched.rda")
tmp2 <- getTime(B_switched,c(3,2),whichValue=c("number.drivers","number.superdrivers"))
save(list=c("tmp","tmp2"),file="tmp.rdata")
