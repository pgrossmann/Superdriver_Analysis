#---------------------------------
### read in simulation results ###
#---------------------------------

n_toCompare <- 10

### reading examples ###

testDir <- file.path("..", "C", "out", "no_intermediate", "test_simus")
testSimuDirs <- list.files(testDir)

# stop("bla")

print("reading in final pops..")

pops <- readPopOfDirs(sapply(testSimuDirs[1:3], function(x) file.path(testDir, x)), "r.+.pop")

# stop("testing")

### partitioning examples ###

print("partitioning..")

pops_up <- getParts2(pops, "upper", n_toCompare)
pops_middle <- getParts2(pops, "middle", n_toCompare)
pops_lower <- getParts2(pops, "lower", n_toCompare)
pops_up_means <- getMeansOf2(pops_up)
pops_up_sds <- getSDsOf2(pops_up)
pops_middle_sds <- getSDsOf2(pops_middle)
pops_middle_means <- getMeansOf2(pops_middle)
pops_lower_means <- getMeansOf2(pops_lower)


print("reading in sub-pops..")

pops_intermediate <- readPopOfDirs(sapply(testSimuDirs[1:2], function(x) file.path(testDir, x)), "t.+.pop")

print("separating sub-pops..")

pops_intermediate_separated <- separateIntermediatePops2(pops_intermediate, 'g')
pops_intermediate_separated2 <- separateIntermediatePops(pops_intermediate[[1]], 'g')
pops_intermediate_separated_up <- getParts3(pops_intermediate_separated, "upper", n_toCompare)
pops_intermediate_separated_up_means <- getMeansOf3(pops_intermediate_separated_up)

pops_intermediate_separated_up_maxs <- getOneValueStatForeach_max3(pops_intermediate_separated_up)
pops_intermediate_separated_up_maxsOfMutations <- getFromOneValueStat3(pops_intermediate_separated_up_maxs, "number.mutations")

pops_intermediate_separated_maxs <- getOneValueStatForeach_max3(pops_intermediate_separated)
pops_intermediate_separated_maxsOfMutations <- getFromOneValueStat3(pops_intermediate_separated_maxs, "number.mutations")

print("get waiting time")
waitingTimes_cancer <- getWaitingTimes_toCancer(pops_intermediate_separated, 2)
waitingTimes_cancer2 <- getWaitingTimes_toCancer_store(pops_intermediate_separated, 2, store=TRUE, cores=1)
waitingTimes_cancer3 <- getWaitingTimes_toCancer(switchPops3(pops_intermediate_separated), 2)
waitingTimes_cancer4 <- getWaitingTime_toCancer(switchPops2(pops_intermediate_separated[[1]]), 2)
waitingTimes_cancer5 <- getWaitingTime_toCancer(pops_intermediate_separated[[1]], c(1,1))
waitingTimes_cancer6 <- getWaitingTime_toCancer(switchPops2(pops_intermediate_separated[[1]]), c(1,1))
waitingTimes_mutator <- getWaitingTimes_toMutator(pops_intermediate_separated)
waitingTimes_mutator2 <- getWaitingTimes_toMutator(pops_intermediate_separated, 2)

print("get one 4mut object")
mut4_test_notWeighted <- get4MutationClasses(pops_intermediate_separated[[1]][[1]][[1]][1:10,], weighted=FALSE)
mut4_test_weighted <- get4MutationClasses(pops_intermediate_separated[[1]][[1]][[1]][1:10,], weighted=TRUE)
mut4_test_weighted2 <- get4MutationClasses3(pops_intermediate_separated, weighted=TRUE)
mut4_test_weighted3 <- get4MutationClasses3_store(pops_intermediate_separated, weighted=FALSE, store=TRUE, cores=1)
mut4_test_weighted4 <- get4MutationClasses_mean2(pops_intermediate_separated[[1]], weighted=TRUE)

print("get N 4mut classes")
Nmut4_test_weighted4 <- getNumber4mutClasses1(mut4_test_weighted4)
Nmut4_test_weighted3 <- getNumber4mutClasses3(mut4_test_weighted3)

print("get probs")
probs_old <- lapply(getProb3(pops_intermediate_separated, "number.mutations", 2), unlist)
probs <- getWaitingTimesProb_toCancer3(pops_intermediate_separated, 4)
probs2 <- getWaitingTimesProb_toCancer3_store(pops_intermediate_separated, 4, store=TRUE, cores=1)

print("analyseMore")
sourcef("analyseMore_functions.R")
waitingTimes_cancer33 <- waitingTimes_cancer3
waitingTimes_cancer33[[2]] <- waitingTimes_cancer33[[1]]
names(waitingTimes_cancer33) <- sapply(names(waitingTimes_cancer33), function(x) sprintf("%sratenormal0.0001",x))
# plotWaitingTimes_twoFixed(waitingTimes_cancer33, mutator=1000, driver=0.005, factor=".*", ratenormal="")

plotFun_varySome_test2(waitingTimes_cancer33, c("mutator","ratenormal"), c("driver"), fix=c("factor"=1.75),
                       ylab="test1", funCenter=function(x) mean(x$repGens), funSd=function(x) sd(x$repGens),
                       epsil=0.0001)

plotFun_varySome_test2(waitingTimes_cancer33, c("mutator","ratenormal"), c("factor"), fix=c("driver"=0.005),
                       ylab="test1", funCenter=function(x) mean(x$repGens), funSd=function(x) sd(x$repGens),
                       epsil=0.01)

plot_varySome.waitingTime(waitingTimes_cancer33, c("mutator","ratenormal"), c("factor"), fix=c("driver"=0.005),
                          ylab="test2", funCenter=function(x) mean(x$repGens), funSd=function(x) sd(x$repGens),
                          epsil=0.01)

plot_varySome.waitingTime.cancer_forReps(waitingTimes_cancer33, c("mutator","ratenormal"), c("factor"), 
                                         fix=c("driver"=0.005), epsil=0.01, k=2)

plot_varySome.waitingTime.cancer_forReps(waitingTimes_cancer33, c("mutator","ratenormal"), c("driver"), 
                                         fix=c("factor"=1.75), epsil=0.01, k=2)


mut4_test_weighted5 <- mut4_test_weighted3
mut4_test_weighted5[[1]] <- mut4_test_weighted3[[1]][[1]]
mut4_test_weighted5[[2]] <- mut4_test_weighted3[[1]][[2]]
names(mut4_test_weighted5) <- sapply(names(mut4_test_weighted5), function(x) sprintf("%sratenormal0.0001",x))

pdf("/tmp/tmp.pdf")
plot_varySome.N4mutClasses.finalPops(getNumber4mutClasses2(mut4_test_weighted5), c("mutator","ratenormal"), 
                                     c("factor"), fix=c("driver"=0.005), epsil=0.001)

plot_varySome.N4mutClasses.finalPops(getNumber4mutClasses2(mut4_test_weighted5), c("mutator","factor"), 
                                     c("driver"), fix=c("ratenormal"=0.0001), epsil=0.0001)

plot_varySome.N4mutClasses.finalPops(getNumber4mutClasses2(mut4_test_weighted5), c("factor","mutator"), 
                                     c("driver"), fix=c("ratenormal"=0.0001), epsil=0.0001)
dev.off()
stop("stop")
plotWaitingTimes_twoFixed_forAllCase("A.+t_mutator.rda")
plotNclasses_twoFixed_forAllCase("C.+4mut.+.rda")

stop("stop before plotting examples")

### plotting examples ###

print("plotting examples")

pdf(toPlots("examples_compareMeans.pdf"))
plotFrequencies(pops[[1]][[1]])

tmp_compareMeans <- list(Up=pops_up_means, Mid=pops_middle_means)
tmp_compareSDs <- list(Up=pops_up_sds, Mid=pops_middle_sds)

scatterPlot(tmp_compareMeans, columns=c("fitness","count"), lLTabs_sd=tmp_compareSDs)
scatterPlot(tmp_compareMeans, columns=c("fitness"))
scatterPlot(tmp_compareMeans, columns=c("count"), histType="highLow")
scatterPlot(tmp_compareMeans, columns=c("fitness"), histType="lowHigh")
scatterPlot(tmp_compareMeans, columns=c("fitness"), histType="lowHigh", lLTabs_sd=tmp_compareSDs)

dev.off()

pdf(toPlots("examples_compareReplicates.pdf"))

tmp_compareReplicates <- list("12"=pops_up[1:2],"3"=pops_up[3])

scatterPlot2(tmp_compareReplicates, columns=c("fitness","count"), levels=3)
scatterPlot2(tmp_compareReplicates, columns=c("fitness"), levels=3)
scatterPlot2(tmp_compareReplicates, columns=c("count"), histType="highLow", levels=3)
scatterPlot2(tmp_compareReplicates, columns=c("fitness"), histType="lowHigh", levels=3)

dev.off()

## try out heatmap
tmpPopsMeans <- getMeansOf2(getParts2(pops, "upper", 100))[[1]][,-1]
require(gplots)
ab <- range(tmpPopsMeans)
bk <- seq(ab[1],ab[2],by=0.1)
cl <- colorpanel(length(bk)-1,"yellow","slateblue4")
pdf(toPlots("example_heatmap.pdf"),width=25,height=30)
heatmap.2(tmpPopsMeans,dendrogram="column",Rowv=F,breaks=bk,col=cl,key=F)
dev.off()

pdf(toPlots("fake_generations.pdf"))
plot(exp(seq(4,9,length=50))*rexp(1:100), col="green", type="l", xlab="Generation", ylab="arbitrary value")
lines(exp(seq(1,10,length=100))*rexp(1:100), col="blue")
dev.off()
