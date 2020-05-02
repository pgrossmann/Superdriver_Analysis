###########################
### 
### performs pre-analysis 
###
##########################

sourcef("load_analyse_store2_functions.R")

## just a convenient way to run stuff independently by uncommenting when needed
execWaitingTime <- FALSE
exec4mutClasses_final <- FALSE
execoneClones_final <- FALSE
exec4mutClasses_inter <- FALSE
execoneClones_inter <- FALSE

execWaitingTime <- TRUE
# exec4mutClasses_final <- TRUE
# execoneClones_final <- TRUE
# exec4mutClasses_inter <- TRUE
# execoneClones_inter <- TRUE

# execWaitingTime_only <- TRUE
# execN4mutClasses_only <- TRUE
# execoneClones_only <- TRUE
# execN4mutClasses_only <- TRUE
# execoneClones_only <- TRUE

# execWaitingTime_only <- FALSE
# execN4mutClasses_only <- FALSE
# execoneClones_only <- FALSE

##########################
### waiting time analysis
##########################

if (execWaitingTime) {
  # waitingTimesA <- getAnalysis.waitingTimes("^A_mutator1000driver0.001factor1", parallel=FALSE)
  #
  #getAnalysis.waitingTimes("^C_mutator1.*factor1.5ratenormal0.00000001")
  
  #stop("special cases done")
  
  ## u=1e-8
  for (i in 1:10) {
    for (j in 1:6) {
      getAnalysis.waitingTimes_KrKs("^B_.*ratenormal0.00000001",kCancer=c(kr=i,ks=j))
      # getAnalysis.waitingTimes_KrKp("^A_.*ratenormal0.00000001",kCancer=c(kp=i,kr=j))
    }
  }
  ## u=1e-7
  for (i in 1:10) {
    for (j in 1:6) {
      # getAnalysis.waitingTimes_KrKs("^B_.*ratenormal0.0000001",kCancer=c(kr=i,ks=j))
      # getAnalysis.waitingTimes_KrKp("^A_.*ratenormal0.0000001",kCancer=c(kp=i,kr=j))
    }
  }
  
#   stop("waiting time special cases done")
  
  ## for u=1e-8
  # waitingTimesA10eminus8 <- 
  #   getAnalysis.waitingTimes("^A_.*ratenormal0.00000001") #mutator1000driver0.1factor1ratenormal0.0000001
  # getAnalysis.waitingTimes("^B_.*ratenormal0.00000001")
  # getAnalysis.waitingTimes("^C_.*ratenormal0.00000001")
  # getAnalysis.waitingTimes("^D_.*ratenormal0.00000001")
  # 
  # ## for u=10e-7
  # waitingTimesA10eminus7 <- 
  #   getAnalysis.waitingTimes("^A_.*ratenormal0.0000001") #mutator1000driver0.1factor1ratenormal0.0000001
  # getAnalysis.waitingTimes("^B_.*ratenormal0.0000001")
  # getAnalysis.waitingTimes("^C_.*ratenormal0.0000001")
  # getAnalysis.waitingTimes("^D_.*ratenormal0.0000001")
  print(proc.time())
  # 
#   if (execWaitingTime_only)
    stop("waiting time done")
}

##############
### Diversity
##############

#---------------------
### final populations 
#---------------------

if (exec4mutClasses_final) {
#   classesD_finalPops <- getAnalysis.4mutClasses_finalPops("^finalPops_D_")
#   classesMeansD_finalPops <- getAnalysis.4mutClassesMeans_finalPops("^finalPops_D_")
  
#   stop("case D finished")
  
  classesA_finalPops <- getAnalysis.4mutClasses_finalPops("^finalPops_A_")
  classesB_finalPops <- getAnalysis.4mutClasses_finalPops("^finalPops_B_")
  classesC_finalPops <- getAnalysis.4mutClasses_finalPops("^finalPops_C_")
  classesD_finalPops <- getAnalysis.4mutClasses_finalPops("^finalPops_D_")
  
  classesMeansA_finalPops <- getAnalysis.4mutClassesMeans_finalPops("^finalPops_A_")
  classesMeansB_finalPops <- getAnalysis.4mutClassesMeans_finalPops("^finalPops_B_")
  classesMeansC_finalPops <- getAnalysis.4mutClassesMeans_finalPops("^finalPops_C_")
  classesMeansD_finalPops <- getAnalysis.4mutClassesMeans_finalPops("^finalPops_D_")
  
  # classesA_finalPops_test <- getAnalysis.4mutClasses_finalPops("finalPops_A_mutator1000driver0.001factor1", parallel=FALSE)
  # classesMeansA_finalPops_test <- getAnalysis.4mutClassesMeans_finalPops("finalPops_A_mutator1000driver0.001factor1", parallel=FALSE)
  
  stop("4mutclasses finalPops")
}

if (execoneClones_final) {
  # oneClonesA_finalPops_test <- getAnalysis.oneClones_finalPops("finalPops_A_mutator1000driver0.001factor1", parallel=FALSE)
  
  ## oneClones on final
  
  oneClonesA_finalPops <- getAnalysis.oneClones_finalPops("^finalPops_A_")
  getAnalysis.oneClones_finalPops("^finalPops_B_")
  getAnalysis.oneClones_finalPops("^finalPops_C_")
  getAnalysis.oneClones_finalPops("^finalPops_D_")
  
  stop("got oneClones finalPops")
}

#----------------------------
### intermediate populations 
#----------------------------

if (exec4mutClasses_inter) {
  goParallel <- TRUE
  # classesA_test <- getAnalysis.4mutClasses("A_mutator1000driver0.001factor1", parallel=FALSE)
  # stop("test")
#   getAnalysis.4mutClassesMeans("^[AB]_.*mutator1000.*driver0.01factor1.5ratenormal0.00000001.*", whichWeight="count")
#   getAnalysis.4mutClassesMeans("^[AB]_.*mutator1000.*driver0.001factor1.5ratenormal0.00000001.*", whichWeight="count")
#   getAnalysis.4mutClassesMeans("^[AB]_.*mutator1000.*driver0.001factor1.5ratenormal0.0000001.*", whichWeight="count")
#   stop("4mut classes_inter special cases done")
#   print("exec4mutClasses_inter")
#   # 1e-8
#   getAnalysis.4mutClassesMeans("^B_.*driver0.01factor2ratenormal.*", whichWeight="count")
#   getAnalysis.4mutClassesMeans("^[AB]_.*mutator1000.*driver0.01factor3ratenormal0.00000001.*", whichWeight="count")
#   getAnalysis.4mutClassesMeans("^[AB]_.*mutator1000.*driver0.01factor1ratenormal0.00000001.*", whichWeight="count")
#   getAnalysis.4mutClassesMeans("^[AB]_.*mutator1000.*driver0.01factor1.5ratenormal0.00000001.*", whichWeight="count")
#   getAnalysis.4mutClassesMeans("^[AB]_.*mutator1000.*driver0.001factor3ratenormal0.00000001.*", whichWeight="count")
#   getAnalysis.4mutClassesMeans("^[AB]_.*mutator1000.*driver0.1factor3ratenormal0.00000001.*", whichWeight="count")
#   getAnalysis.4mutClassesMeans("^[AB]_.*mutator1000.*driver0.1factor1ratenormal0.00000001.*", whichWeight="count")
#   getAnalysis.4mutClassesMeans("^[AB]_.*mutator1000.*driver0.001factor1ratenormal0.00000001.*", whichWeight="count")
#   getAnalysis.4mutClassesMeans("^[AB]_.*mutator1000.*driver0.001factor1.5ratenormal0.00000001.*", whichWeight="count")
#   # 1e-7
#   getAnalysis.4mutClassesMeans("^[AB]_.*mutator1000.*driver0.01factor3ratenormal0.0000001.*", whichWeight="count")
#   getAnalysis.4mutClassesMeans("^[AB]_.*mutator1000.*driver0.01factor1ratenormal0.0000001.*", whichWeight="count")
#   getAnalysis.4mutClassesMeans("^[AB]_.*mutator1000.*driver0.01factor1.5ratenormal0.0000001.*", whichWeight="count")
#   getAnalysis.4mutClassesMeans("^[AB]_.*mutator1000.*driver0.001factor3ratenormal0.0000001.*", whichWeight="count")
#   getAnalysis.4mutClassesMeans("^[AB]_.*mutator1000.*driver0.1factor3ratenormal0.0000001.*", whichWeight="count")
#   getAnalysis.4mutClassesMeans("^[AB]_.*mutator1000.*driver0.1factor1ratenormal0.0000001.*", whichWeight="count")
#   getAnalysis.4mutClassesMeans("^[AB]_.*mutator1000.*driver0.001factor1ratenormal0.0000001.*", whichWeight="count") 
#   getAnalysis.4mutClassesMeans("^[AB]_.*mutator1000.*driver0.001factor1.5ratenormal0.0000001.*", whichWeight="count")
#   stop("B 4mutClassesMeans special cases done")

  getAnalysis.4mutClasses("^[AB]_.*mutator1000.*driver0.01factor2.6ratenormal0.00000001.*", parallel=F) ## EDIT August, 2, 2017
  # getAnalysis.4mutClassesMeans("^[AB]_.*mutator1000.*driver.*factor2ratenormal0.00000001.*", whichWeight="count")

stop("new special cases done")
  
  getAnalysis.4mutClassesMeans("^A_", whichWeight="count")
  getAnalysis.4mutClassesMeans("^B_", whichWeight="count")
  getAnalysis.4mutClassesMeans("^C_", whichWeight="count")
  getAnalysis.4mutClassesMeans("^D_", whichWeight="count")
  
  getAnalysis.4mutClasses("^A_")
  getAnalysis.4mutClasses("^B_")
  getAnalysis.4mutClasses("^C_")
  getAnalysis.4mutClasses("^D_")
  
  # oneClonesA_test <- getAnalysis.oneClones("A_mutator1000driver0.001factor1", parallel=FALSE)
  stop("4mut classes inter")
}

if (execoneClones_inter) {
  getAnalysis.oneClones("^A_")
  getAnalysis.oneClones("^B_")
  getAnalysis.oneClones("^C_")
  getAnalysis.oneClones("^D_")
  
  stop("oneclones inter")
}

# stop("overall analysis finished")
