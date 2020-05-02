sourcef("analyseMore_functions.R")

##################################
# auxillary functions

tmp_getObjects <- function(specificCase, parallel, cores) {
  files <- grep(specificCase, list.files(intermediateDir), value=T) 

  fileNames <- sub("(.+)\\.rda", "\\1", files)
  if (!parallel) {
    objects <- lapply(fileNames, function(fileName) {
      print("before loadIfNot in tmp_getObjects()")
      print(fileName)
      loadIfNot(file.path(intermediateDir,sprintf("%s.rda",fileName)))
    })
  } else {
    require(parallel)
    objects <- mclapply(fileNames, function(fileName) {
      print("before loadIfNot in tmp_getObjects()")
      print(fileName)
      loadIfNot(file.path(intermediateDir,sprintf("%s.rda",fileName)))
    }, mc.cores=cores)
  }
  names(objects) <- fileNames
  return(objects)
}

debsub <- function(x) deparse(substitute(x))
gett <- function(x) if (!is.null(x)) get(x) else x
##########################
### waiting time analysis
##########################

plotWaitingTimesResults <- function(case, onPlot, onX, fix, kCancer, epsil=0.01, parallel=goParallel, 
                                    cores=nCores, xloga=NULL, yloga=NULL) {
#   # for all waiting time res objects of case
#   files <- grep(sprintf("%s.*waitingTimeRes.*\\.rda",case), list.files(intermediateDir), value=T) 
# #   files <- sapply(files, function(x) file.path(intermediateDir,x))
#   fileNames <- sub("(.+)\\.rda", "\\1", files)
#   objects <- lapply(fileNames, function(fileName) {
#     loadIfNot(file.path(intermediateDir,sprintf("%s.rda",fileName)))
#     return(get(fileName))
#   })
#   names(objects) <- fileNames
  
  objects <- tmp_getObjects(sprintf("%s.*waitingTimeResK%s.*\\.rda",case,kCancer), parallel=parallel, cores=cores)
  
  print("got waiting time files")
#   print(fileNames)

  if (length(objects) == 0) {
    print(sprintf("length objects 0, with k=%s. nothing to plot!",kCancer))
    return(NULL)
  }
  
  if (is.null(xloga))
    pdfName <- sprintf("%swaitingTimeResK%s.pdf",case,kCancer)
  else
    pdfName <- sprintf("%swaitingTimeResK%s_x%s.pdf",case,kCancer,as.character(substitute(xloga)))
  
#   pdf(file.path(plotDir,pdfName))
  objects_cancerWaitingTime <- lapply(objects, function(x) x$cancerWaitingTimeSwitched)
  objects_mutatorWaitingTime <- lapply(objects, function(x) x$mutatorWaitingTimeSwitched)
#   xloga <- gett(xloga)
#   yloga <- gett(yloga)
  plot_varySome.waitingTime.cancer_forReps(objects_cancerWaitingTime, onPlot=onPlot, onX=onX, 
                                           fix=fix, epsil=epsil, k=kCancer, xloga=xloga, yloga=yloga)
  plot_varySome.waitingTime.mutator_forReps(objects_mutatorWaitingTime, onPlot=onPlot, onX=onX,
                                            fix=fix, epsil=epsil)
#   dev.off()
}

## k = 20
# forAllofOneParam_printTo1pdf <- function
plotWaitingTime4diffK2plot <- function(k_c) {
  pdf(file.path(plotDir,sprintf("test_%swaitingTimeResK%s.pdf","A_",k_c)))
  plotWaitingTimesResults("A_", fix=c("mutator"=1000), onPlot=c("factor","ratenormal"), 
                          kCancer=k_c, onX="driver", epsil=0.001)
  dev.off()
  
  pdf(file.path(plotDir,sprintf("test_%swaitingTimeResK%s.pdf","B_",k_c)))
  for (u in c(0.0000001,0.00000001)) {
    plotWaitingTimesResults("B_", fix=c("mutator"=1000,"ratenormal"=u), onPlot=c("factor"),
                            kCancer=k_c, onX="driver", epsil=0.001)
  }
  dev.off()
  
  pdf(file.path(plotDir,sprintf("test_%swaitingTimeResK%s.pdf","C_",k_c)))
  for (m in 5:1) {
    plotWaitingTimesResults("C_", fix=c("ratenormal"=0.00000001,"mutator"=m), onPlot=c("factor"), 
                            kCancer=k_c, onX="driver", epsil=0.001)
    plotWaitingTimesResults("C_", fix=c("ratenormal"=0.00000001,"mutator"=m), onPlot=c("factor"), 
                            kCancer=k_c, onX="driver", epsil=0.001)
  }
  dev.off()
  
  pdf(file.path(plotDir,sprintf("test_%swaitingTimeResK%s.pdf","D_",k_c)))
  for (u in c(0.0000001,0.00000001)) {
    plotWaitingTimesResults("D_", fix=c("ratenormal"=u,"factor"=1), onPlot=c("mutator"), 
                            kCancer=k_c, onX="driver", epsil=0.001)
  }
  dev.off()
}

plotWaitingTime4diffK2plot(10)
plotWaitingTime4diffK2plot(20)
stop("waitingTime plots")
# 
# stop("test")
# 
# # 
# # ## k = 10
# plotWaitingTimesResults("A_", fix=c("mutator"=1000), onPlot=c("factor","ratenormal"), 
#                         kCancer=10, onX="driver", epsil=0.001)
# plotWaitingTimesResults("B_", fix=c("mutator"=1000), onPlot=c("factor","ratenormal"),
#                         kCancer=10, onX="driver", epsil=0.001)
# plotWaitingTimesResults("C_", fix=c("ratenormal"=0.00000001), onPlot=c("mutator","factor"), 
#                         kCancer=10, onX="driver", epsil=0.001, xloga=log10)
# plotWaitingTimesResults("D_", fix=c("ratenormal"=0.00000001), onPlot=c("mutator","factor"), 
#                         kCancer=10, onX="driver", epsil=0.001)

##############
### Diversity
##############

#---------------------
### final populations 
#---------------------

### mutation classes ###

plot4mutClassResults.final <- function(case, onPlot, onX, fix, epsil=0.01, parallel=goParallel, 
                                       cores=nCores, xloga=NULL, yloga=NULL) {
#   # for all waiting time res objects of case
#   files <- grep(sprintf("%s.*4mutClassesResFinalPops.*\\.rda",case), list.files(intermediateDir), value=T) 
#   #   files <- sapply(files, function(x) file.path(intermediateDir,x))
#   files <- files
#   fileNames <- sub("(.+)\\.rda", "\\1", files)
# #   
#   if (!parallel) {
#     objects <- lapply(fileNames, function(fileName) {
#           print("before loadIfNot")
#           print(fileName)
#           loadIfNot(file.path(intermediateDir,sprintf("%s.rda",fileName)))
#       #     print("lapplying")
#       #     return(get(fileName))
#         })
#   } else {
#     require(parallel)
#     objects <- mclapply(fileNames, function(fileName) {
#       print("before loadIfNot")
#       print(fileName)
#       loadIfNot(file.path(intermediateDir,sprintf("%s.rda",fileName)))
#       #     print("lapplying")
#       #     return(get(fileName))
#     }, mc.cores=cores)
#   }
#   names(objects) <- fileNames
  
  objects <- tmp_getObjects(sprintf("%s.*4mutClassesResFinalPops.*\\.rda",case), parallel=parallel, cores=cores)
  
  if (length(objects) == 0) {
    print(sprintf("length objects 0, with 4mutClasses of case %s. nothing to plot!",case))
    return(NULL)
  }
  
  print("got N classes files")
#   print(fileNames)
  if (is.null(xloga))
    pdfName <- sprintf("%s4mutClassesResFinalPops.pdf",case)
  else
    pdfName <- sprintf("%s4mutClassesResFinalPops_xlog%s.pdf",case,as.character(substitute(xloga)))
  
#   pdf(file.path(plotDir,pdfName))
#   objects <- lapply(objects, function(x) x$fourMutClasses)
  objects_numberClasses <- lapply(objects, function(x) x$N4mutClasses)
#   objects_numberClasses <- getNumber4mutClasses2(objects)
#   print(objects_numberClasses)
  print("about to plot")
  plot_varySome.N4mutClasses.finalPops(objects_numberClasses, onPlot=onPlot, onX=onX, fix=fix,
                                       epsil=epsil, plotSd=T, xloga=xloga, yloga=yloga)
#   dev.off()
  objects_numberClasses
}


# finalPops_A_mutator1000driver0.001factor1ratenormal0.00000001_4mutClassesResFinalPops
# finalPops_A_mutator1000driver0.001factor1ratenormal0.00000001_4mutClassesResFinalPops
# goParallel=FALSE
plot4mutClassResults.final("A_", fix=c("ratenormal"=0.00000001), onPlot=c("factor","mutator"), 
                                   onX="driver", epsil=0.001)
plot4mutClassResults.final("B_", fix=c("ratenormal"=0.00000001), onPlot=c("factor","mutator"), 
                                   onX="driver", epsil=0.001)
plot4mutClassResults.final("B_", fix=c("ratenormal"=0.00000001), onPlot=c("factor","mutator"), 
                           onX="driver", epsil=0.001, xloga=log)
plot4mutClassResults.final("C_", fix=c("ratenormal"=0.00000001), onPlot=c("mutator","factor"), 
                                  onX="driver", epsil=0.001)
plot4mutClassResults.final("D_", fix=c("ratenormal"=0.00000001), onPlot=c("mutator","factor"), 
                                   onX="driver", epsil=0.001)

### average number of (super)drivers ###

plotOneClones_highestFreq.final <- function(case, onPlot, onX, fix, epsil=0.01, parallel=goParallel, 
                                            cores=nCores, xloga=NULL, yloga=NULL) {  
  objects <- tmp_getObjects(sprintf("%s.*oneClonesFinalPopsRes.*\\.rda",case), parallel=parallel, cores=cores)
  
#   print("the oneConesFinalRes files")
#   print(fileNames)
  
  if (length(objects) == 0) {
    print(sprintf("length objects 0, with oneClones of case %s. nothing to plot!",case))
    return(NULL)
  }
  
  #   objects <- lapply(objects, function(x) x$fourMutClasses)
  objects_highest <- lapply(objects,  "[[", "highest") # results in list of replicates of highest
  # object_highest is list (of simus) of list (of replicates) of vectors (leading clone)
  objects_highest_superdrivers <- lapply(objects_highest, function(rep) sapply(rep, "[[", "number.superdrivers"))
  # results in list of simus of vectors of (super)driver numbers for each replicate
  objects_highest_drivers <- lapply(objects_highest, function(rep) sapply(rep, "[[", "number.drivers"))
  objects_highest_mutators <- lapply(objects_highest, function(rep) sapply(rep, "[[", "number.mutators"))
  #   objects_numberClasses <- getNumber4mutClasses2(objects)
  #   print(objects_numberClasses)
  
  ## start pdf here ##
  if (is.null(xloga))
    pdfName <- sprintf("%soneClonesFinalPopsRes_leading.pdf",case)
  else
    pdfName <- sprintf("%soneClonesFinalPopsRes_leading_xlog%s.pdf",case,as.character(substitute(xloga)))
  
#   pdf(file.path(plotDir,pdfName))  
  print("about to plot")
#   xloga <- gett(xloga)
#   yloga <- gett(yloga)
  plot_varySome.oneClonesMutations.finalPops(objects_highest_superdrivers, 
                                             onPlot=onPlot, onX=onX, fix=fix, 
                                             whichMutation="superdrivers", 
                                             epsil=epsil, plotSd=T,
                                             xloga=xloga, yloga=yloga)
  plot_varySome.oneClonesMutations.finalPops(objects_highest_drivers, 
                                             onPlot=onPlot, onX=onX, fix=fix, 
                                             whichMutation="drivers", 
                                             epsil=epsil, plotSd=T,
                                             xloga=xloga, yloga=yloga)
  plot_varySome.oneClonesMutations.finalPops(objects_highest_mutators, 
                                             onPlot=onPlot, onX=onX, fix=fix, 
                                             whichMutation="mutators", 
                                             epsil=epsil, plotSd=T,
                                             xloga=xloga, yloga=yloga)
#   dev.off()
  objects_highest
}

tmp11 <- plotOneClones_highestFreq.final("A_", fix=c("ratenormal"=0.00000001), onPlot=c("factor","mutator"), 
                                   onX="driver", epsil=0.001)
tmp11 <- plotOneClones_highestFreq.final("A_", fix=c("ratenormal"=0.00000001), onPlot=c("factor","mutator"), 
                                         onX="driver", epsil=0.001)
tmp22 <- plotOneClones_highestFreq.final("B_", fix=c("ratenormal"=0.00000001), onPlot=c("factor","mutator"), 
                                   onX="driver", epsil=0.001)
tmp33 <- plotOneClones_highestFreq.final("C_", fix=c("ratenormal"=0.00000001), onPlot=c("mutator","factor"), 
                                   onX="driver", epsil=0.001)
tmp44 <- plotOneClones_highestFreq.final("D_", fix=c("ratenormal"=0.00000001), onPlot=c("mutator","factor"), 
                                   onX="driver", epsil=0.001)

#----------------------------
### intermediate populations 
#----------------------------
