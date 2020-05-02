sourcef("analyseMore_functions.R")
sourcef("modelFitting.R") # in plotList

##################################
# auxillary functions

tmp_getObjects <- function(specificCase, parallel, cores, thedir=intermediateDir) {
  files <- grep(specificCase, list.files(thedir), value=T) 
  
  fileNames <- sub("(.+)\\.rda", "\\1", files)
  if (!parallel) {
    objects <- lapply(fileNames, function(fileName) {
      print("before loadIfNot in tmp_getObjects()")
      print(fileName)
      loadIfNot(file.path(thedir,sprintf("%s.rda",fileName)))
    })
  } else {
    require(parallel)
    objects <- mclapply(fileNames, function(fileName) {
      print("before loadIfNot in tmp_getObjects()")
      print(fileName)
      loadIfNot(file.path(thedir,sprintf("%s.rda",fileName)))
    }, mc.cores=cores)
  }
  names(objects) <- fileNames
  return(objects)
}

debsub <- function(x) deparse(substitute(x))
gett <- function(x) if (!is.null(x)) get(x) else x

plotFun4thisParamValues <- function(cases, casesPlotnames, plotf, ...) {
  casesPlotnames <- sapply(casesPlotnames, function(x) sprintf("%s.pdf",x))
  
  if ("A"%in%names(cases)) {
    pdf(file.path(plotDir,casesPlotnames["A"]))
    plotf(cases["A"], fix=c("mutator"=1000,"factor"=1), onPlot=c("ratenormal"), 
          onX="driver", epsil=0.001, plotApproximate=FALSE, ...)
    
    dev.off()
  } 
  
  if ("B"%in%names(cases)) {
    pdf(file.path(plotDir,casesPlotnames["B"]))
    for (u in c(0.0000001,0.00000001)) {
      plotf(cases["B"], fix=c("mutator"=1000,"ratenormal"=u), onPlot=c("factor"),
            onX="driver", epsil=0.001, plotApproximate=FALSE, ...)
    }
    dev.off()
  }  
    
  if ("C"%in%names(cases)) {
    pdf(file.path(plotDir,casesPlotnames["C"]))
    for (m in 5:1) {
      plotf(cases["C"], fix=c("ratenormal"=0.00000001,"mutator"=m), onPlot=c("factor"), 
            onX="driver", epsil=0.001, plotApproximate=FALSE, ...)
      plotf(cases["C"], fix=c("ratenormal"=0.00000001,"mutator"=m), onPlot=c("factor"), 
            onX="driver", epsil=0.001, plotApproximate=FALSE,...)
    }
    dev.off()
  }
  
  if ("D"%in%names(cases)) {
    pdf(file.path(plotDir,casesPlotnames["D"]))
    for (u in c(0.0000001,0.00000001)) {
      plotf(cases["D"], fix=c("ratenormal"=u,"factor"=1), onPlot=c("mutator"), 
            onX="driver", epsil=0.001, plotApproximate=FALSE, ...)
    }
    dev.off()
  }
   
}

##########################
### waiting time analysis
##########################

plotWaitingTimesResults <- function(case, onPlot, onX, fix, kCancer, epsil=0.01, parallel=goParallel, 
                                    cores=nCores, xloga=NULL, yloga=NULL, ...) {
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
  approx2 <- grepl("A_",pdfName)
  plot_varySome.waitingTime.cancer_forReps(objects_cancerWaitingTime, onPlot=onPlot, onX=onX, 
                                           fix=fix, epsil=epsil, k=kCancer, xloga=xloga, 
                                           yloga=yloga, approx2=approx2, kCancer=kCancer, ...)
#   if (grepl("A_",pdfName)) {
#     xx <- c(0.001,0.005,0.01,0.05,0.1)
#     yy <- tk_niko(k=kCancer,s_r=xx,N_f=10^9,N_0=10^6,u=1e-8,d_r=100,d_s=10)
#     plot(xx,yy,col="red",type="l",xaxt="n")
#     yy <- tk_niko(k=kCancer,s_r=xx,N_f=10^9,N_0=10^6,u=1e-7,d_r=100,d_s=10)
#     lines(xx,yy,col="black")
#     ##
#     xx <- c(0.001,0.005,0.01,0.05,0.1)
#     yy <- tk_new(k=kCancer,s_r=xx,N_f=10^9,N_0=10^6,u=1e-8,d_r=100,d_s=10)
#     lines(xx,yy,col="red",lty=2)
#     yy <- tk_new(k=kCancer,s_r=xx,N_f=10^9,N_0=10^6,u=1e-7,d_r=100,d_s=10)
#     lines(xx,yy,col="black",lty=2)
#     axis(1, at=xx)
#     legend("topright",legend=c("tk_niko 1e-8","tk_niko 1e-7","tk_new 1e-8","tk_new 1e-7"),
#            col=c("black","red","black","red"), lty=c(1,1,2,2))
#   }
  plot_varySome.waitingTime.mutator_forReps(objects_mutatorWaitingTime, onPlot=onPlot, onX=onX,
                                            fix=fix, epsil=epsil, ...)
  #   dev.off()
}

## k = 20
# forAllofOneParam_printTo1pdf <- function
plotWaitingTime4diffK2plot <- function(k_c) {
  cases <- sapply(LETTERS[1:4], function(x) paste0(x,"_"))
  casesPlotnames <- sapply(cases, function(x) sprintf("%swaitingTimeResK%s",x,k_c))
  tt <- plotFun4thisParamValues(cases=cases, casesPlotnames=casesPlotnames, 
                                plotWaitingTimesResults, kCancer=k_c)
}

##############################
### plot the waiting times ###
# plotWaitingTime4diffK2plot(10)
# plotWaitingTime4diffK2plot(20)
##############################

##############################
### 2d waiting time ###

plotWaitingTimesResults_2d <- function(case, onPlot, onX, fix, kCancer, epsil=0.01, parallel=goParallel, 
                                    cores=nCores, xloga=NULL, yloga=NULL, approxB=TRUE, ...) {
  
  
  objects <- tmp_getObjects(sprintf("%s.*waitingTimeResKr%sKs%s.*\\.rda",case,kCancer["kr"],kCancer["ks"]), 
                            parallel=parallel, cores=cores)
  
  print("got waiting time files")
  #   print(fileNames)
  
  if (length(objects) == 0) {
    print(sprintf("length objects 0, with k=%s. nothing to plot!",kCancer))
    return(NULL)
  }
  
#   if (is.null(xloga))
#     pdfName <- sprintf("%swaitingTimeResK%s.pdf",case,kCancer)
#   else
#     pdfName <- sprintf("%swaitingTimeResK%s_x%s.pdf",case,kCancer,as.character(substitute(xloga)))
  
  #   pdf(file.path(plotDir,pdfName))
  objects_cancerWaitingTime <- lapply(objects, function(x) x$cancerWaitingTimeSwitched)
  plot_varySome.waitingTime_2d.cancer_forReps(objects_cancerWaitingTime, onPlot=onPlot, onX=onX, 
                                           fix=fix, epsil=epsil, xloga=xloga, plotLegend.topRight=T,
                                           yloga=yloga, approxB=approxB, kCancer=kCancer, ...)
  
}

#####################################################
########### EXECUTE WAITING TIME ANALYSIS ###########  

### The process is as follows:
### First, execute this pipeline INCLUDING the model building with
### errorCorrection.lm=FALSE (DONT FORGET TO ADAPT THE FILENAME PDF),
### Second, re-run (without the model building ideally) with errorCorrection.lm=TRUE

for (kr in 1:10) {
  for (ks in 1:6) {
    pdf(file.path(plotDir,sprintf("B_waitingTimeResKr%sKs%s_errorCorrection.lm.pdf",kr,ks)))
    for (u in c(0.00000001)) { #, 0.0000001)) { # EDIT: AUGUST 3, 2018
      aggregatedPlottedWaitingTimes_u <- list()
      res <- plotWaitingTimesResults_2d("B_", fix=c("mutator"=1000,"ratenormal"=u), onPlot=c("factor"),
                                        onX="driver", epsil=0.001, kCancer=c(kr=kr,ks=ks),
                                        errorCorrection.lm=TRUE, yloga=log)
      # EDIT: AUGUST 3, 2018
      save(res, file=file.path(intermediateDir, "plottedValues_waitingTime", sprintf("B_waitingTimeResKr%sKs%s_plottedValues.rda",kr,ks)))
    }
    dev.off()
  }
}

getErrorCorrectionModel <- function() {
  dir <- file.path(intermediateDir, "plottedValues_waitingTime")
  files <- list.files(dir, full.names=T)
  waitingTimes <- lapply(files, function(file) {
    load(file)
    res
  })
  names(waitingTimes) <- basename(files)
  ## calculate the diff
  waitingTimes_diff <- lapply(waitingTimes, function(x) {
    names_theoretical <- names(x$plottedWaitingTimes$theoretical$yy)
    names_simulated <- names(x$plottedWaitingTimes$simulated$yy)
    ins <- intersect(names_theoretical, names_simulated)
    dif <- x$plottedWaitingTimes$simulated$yy[ins] - x$plottedWaitingTimes$theoretical$yy[ins]
    names(dif) <- ins
    dif
  })
  ## Convert list to Matrix
  totalDiffsList <- unlist(waitingTimes_diff)
  Ks.vals <- gsub(".*Ks([0-9]+).*", "\\1", names(totalDiffsList))
  Kr.vals <- gsub(".*Kr([0-9]+).*", "\\1", names(totalDiffsList))
  driver.vals <- gsub(".*driver([0-9\\.]+).*", "\\1", names(totalDiffsList))
  factor.vals <- gsub(".*factor([0-9\\.]+).*", "\\1", names(totalDiffsList))
#   browser()
  df_diffs <- data.frame(Diff=as.vector(totalDiffsList), "Ks"=Ks.vals, "Kr"=Kr.vals, "Driver"=driver.vals, "Factor"=factor.vals)
  df_diffs <- as.data.frame(apply(df_diffs, 2, as.numeric))
  # # EDIT NOVEMBER 2, 2018: HARMONIZE FITNESS LANDSCAPE
  # df_diffs$Factor <- df_diffs$Factor*df_diffs$Driver 
  df_diffs$TotalMutations <- df_diffs$Ks+df_diffs$Kr
  df_diffs$DriverFactorRatio <- df_diffs$Driver/df_diffs$Factor
  df_diffs$FactorDriverRatio <- df_diffs$Factor/df_diffs$Driver
  df_diffs$Superdriver <- df_diffs$Factor*df_diffs$Driver
  
  mod_1 <- lm(Diff ~ ., data=df_diffs[, c("Diff", "Ks"), drop=F])
  mod_2 <- lm(Diff ~ ., data=df_diffs[, c("Diff", "Kr"), drop=F])
  mod_3 <- lm(Diff ~ ., data=df_diffs[, c("Diff", "Ks", "Kr"), drop=F])
  mod_4 <- lm(Diff ~ ., data=df_diffs[, c("Diff", "Driver"), drop=F])
  mod_5 <- lm(Diff ~ ., data=df_diffs[, c("Diff", "Factor"), drop=F])
  mod_6 <- lm(Diff ~ ., data=df_diffs[, c("Diff", "Driver", "Factor"), drop=F])
  mod_7 <- lm(Diff ~ ., data=df_diffs[, c("Diff", "TotalMutations"), drop=F])
  mod_8 <- lm(Diff ~ ., data=df_diffs[, c("Diff", "TotalMutations", "Driver"), drop=F])
  mod_9 <- lm(Diff ~ ., data=df_diffs[, c("Diff", "TotalMutations", "Factor"), drop=F])
  mod_10 <- lm(Diff ~ ., data=df_diffs[, c("Diff", "TotalMutations", "Driver", "Factor"), drop=F])
  mod_11 <- lm(Diff ~ ., data=df_diffs[, c("Diff", "Ks", "Kr", "Driver", "Factor"), drop=F])
  mod_12 <- lm(Diff ~ ., data=df_diffs[, c("Diff", "TotalMutations", "Ks", "Kr", "Driver", "Factor"), drop=F])
  mod_13 <- lm(Diff ~ ., data=df_diffs[, c("Diff", "DriverFactorRatio", "Ks", "Kr", "Driver", "Factor"), drop=F])
  mod_14 <- lm(Diff ~ ., data=df_diffs[, c("Diff", "DriverFactorRatio"), drop=F])
  mod_15 <- lm(Diff ~ ., data=df_diffs[, c("Diff", "FactorDriverRatio", "Ks", "Kr", "Driver", "Factor"), drop=F])
  mod_16 <- lm(Diff ~ ., data=df_diffs[, c("Diff", "FactorDriverRatio"), drop=F])
  mod_17 <- lm(Diff ~ ., data=df_diffs[, c("Diff", "FactorDriverRatio", "Ks", "Kr"), drop=F])
  mod_18 <- lm(Diff ~ ., data=df_diffs[, c("Diff", "DriverFactorRatio", "Ks", "Kr"), drop=F])
  mod_19 <- lm(Diff ~ ., data=df_diffs[, c("Diff", "Factor", "Driver",  "Ks", "Kr"), drop=F])
  mod_20 <- lm(Diff ~ ., data=df_diffs[, c("Diff", "Factor", "Driver"), drop=F])
  mod_21 <- lm(Diff ~ ., data=df_diffs[, c("Diff", "Driver", "Superdriver"), drop=F])

  mods <- list("Ks"=mod_1, "Kr"=mod_2, "KsKr"=mod_3, 
               "Driver"=mod_4, "Factor"=mod_5, "DriverFactor"=mod_6, 
               "TotalMutations"=mod_7, "TotalMutationsDriver"=mod_8, "TotalMutationsFactor"=mod_9, "TotalMutationsDriverFactor"=mod_10,  
               "KsKrDriverFactor"=mod_11, "TotalMutationsKsKrDriverFactor"=mod_12, 
               "DriverFactorRatioKsKrDriverFactor"=mod_13, "DriverFactorRatio"=mod_14,
               "FactorDriverRatioKsKrDriverFactor"=mod_15, "FactorDriverRatio"=mod_16,
               "FactorDriverRatioKsKr"=mod_17, "DriverFactorRatioKsKr"=mod_18,
               "FactorDriverKsKr"=mod_19, "FactorDriver"=mod_20, "DriverSuperdriver"=mod_21)
  mods_adjRs <- lapply(mods, function(mod) summary(mod)$adj.r.squared)
  list(mods=mods, mods_adjRs=mods_adjRs, data=df_diffs)
}

waitingTimeErrorCorrectionModels <- getErrorCorrectionModel()
save(waitingTimeErrorCorrectionModels, file="waiting-time-error-correction-models.rda")

### EDIT SEPTEMBER 2018
### Conclusion:
### Driver per Factor is the highest predictor on non-log-scale
### Factor per Driver is not predictive on non-log scale
### Driver per Factor is not predictive on log-scale
### Factor per Driver is predictive on log-scale
### Highest r-squared with Driver per Factor on log-scale
### 
### Also, the overall r-squared does not improve on log-scale, only the figures look qualitatively better.

stop("tset")

#####################################################

aggregateN4mutClassesMean <- function() {
  n4MutClassesMean_files <- list.files(intermediateDir, pattern="^B_.+N4mutClassesMean\\.rda", full.names=T)
  n4MutClassesMean_objects <- lapply(n4MutClassesMean_files, function(file) {
    name <- gsub("\\.rda", "", basename(file))
    if (! name%in%ls(globalenv())) {
      print(sprintf("Loading %s ...", file))
      load(file, envir=globalenv())
    }
    return(get(name, envir=globalenv()))
  }) 
  
  browser()
}

aggregatePlottedWaitingTimeValues <- function() {
  
}

aggregatePlottedWaitingTimeValues()

aggregateN4mutClassesMean()

##############################################################
# stop("waitingTime plots")
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

#--------------
### heatmap ###
#--------------

#'@param cases character string, e.g. "[ABC]"
computeANDplotHeatmap <- function(cases, heatmap_value_fun, heatmap_text_fun, heatmap_plotFun,
                                  yParam="driver", xParams=c("mutator","factor"), 
                                  fixedParam=c("ratenormal"=0.00000001), k=k_cancer, plotLegend=TRUE, 
                                  zeroCol="lightblue", oneCol="green", naCol="gray", noteCol="red",
                                  Acol="pink", Bcol="blue", Ccol="orange", Dcol="yellow",
                                  uniqueParams=NULL, heatmap_realValue_fun=NULL, heatmap_realText_fun=NULL, ...) {
  tmp1 <- getWaitingTimeTable(cases, k, heatmap_value_fun, heatmap_text_fun, 
                              yParam=yParam, xParams=xParams, fixedParam=fixedParam,
                              uniqueParams=uniqueParams, heatmap_realValue_fun=heatmap_realValue_fun,
                              heatmap_realText_fun=heatmap_realText_fun)  
  if (is.null(tmp1)) # dirty fix
    return(tmp1)
  heatmap_plotFun(tmp1[[1]],tmp1[[2]], margins=c(9,9), cexRow=1.9, cexCol=1.4,
                  zeroCol=zeroCol, oneCol=oneCol, naCol=naCol, noteCol=noteCol, ...)
  legend("top", legend=c("Has reached cancer", "Has not reached cancer",
                             "* Has reached mutator phenotype", "Data missing",
                             "driver", "driver+superdriver", "driver+superdriver+mutator",
                             "driver+mutator"),
         fill=c(oneCol, zeroCol, noteCol, naCol,
                Acol, Bcol, Ccol, Dcol),
         ncol=2)  
  if (!is.null(heatmap_realValue_fun) && !is.null(heatmap_realText_fun)) {
    return(tmp1)
  } else return(NULL)
}

computeANDplotHeatmap.binary <- function(cases, ...) {
  tmp1 <- computeANDplotHeatmap(cases, getHeatmapEntry_binary, getHeatmapOverlayText_binary, plotHeatmap.binary,
                                heatmap_realValue_fun=getHeatmapEntry_smallestGeneration, 
                                heatmap_realText_fun=getHeatmapOverlay_smallestGeneration, ...)
  return(tmp1)
}

computeANDplotHeatmap.binary_2pdf <- function(cases, yParam="driver", xParams=c("mutator","factor"), 
                                              fixedParam=c("ratenormal"=0.00000001), k=k_cancer,
                                              uniqueParams=NULL, ...) {
  options("scipen"=1, "digits"=10)
  caseTitle <-  sprintf("%s_k%s_%s%s.pdf",cases,k,names(fixedParam)[1],fixedParam[1])
  caseTitle <- shortSimulationName(caseTitle)
  options("scipen"=100, "digits"=10)
  pdf(file.path(plotDir,sprintf("heatmap_%s",caseTitle)), width=14, height=9)
  mmain <- sub(".*k(.+)_u(.+)\\.pdf", "k=\\1 u=\\2", caseTitle)
  tmp1 <- computeANDplotHeatmap.binary(cases, yParam=yParam, xParams=xParams, fixedParam=fixedParam,k=k,
                               main=mmain, uniqueParams=uniqueParams, ...)
  dev.off()
  if (!is.null(tmp1)) {
    write.table(tmp1[[3]],file=sprintf("realValues_%s.csv",caseTitle), sep="\t", quote=F,
                row.names=T, col.names=T)
    write.table(tmp1[[4]],file=sprintf("realOverlayText_%s.csv",caseTitle), sep="\t", quote=F, 
                row.names=T, col.names=T)
  }
}

factors <- seq(1,3, by=0.5)
drivers <- c(0.001, 0.005, 0.01, 0.05, 0.1)
mutators <- c(1000,1:5)

uniqueParams <- lapply(list("driver"=drivers, "mutator"=mutators, "factor"=factors), as.character)

#########################################################################################################
### plot heatmaps ###
# computeANDplotHeatmap.binary_2pdf("[ABCD]", k=10, uniqueParams=uniqueParams)
# computeANDplotHeatmap.binary_2pdf("[ABCD]", k=20, uniqueParams=uniqueParams)
# computeANDplotHeatmap.binary_2pdf("[ABCD]", k=10, fixedParam=c("ratenormal"=0.0000001), uniqueParams=uniqueParams)
# computeANDplotHeatmap.binary_2pdf("[ABCD]", k=20, fixedParam=c("ratenormal"=0.0000001), uniqueParams=uniqueParams)
# stop("test")
#########################################################################################################

#---------------------
### mutational waves
#---------------------

#########################################
sourcef("plotPlentyWaves_publication.R")
#########################################

###################

plotMutWaves(pdfPath=file.path(plotDir,"wavesk=10.pdf"), 
             variants=list(drivers=0:18, superdrivers=0:7,
                           mutators=NA, passengers=NA),
             K=10)

stop("HELLLOOOOOOOOOOOOO")
plotMutWaves(pdfPath=file.path(plotDir,"wavesk=20.pdf"), 
             variants=list(drivers=0:18, superdrivers=0:7,
                           mutators=NA, passengers=NA),
             K=20)

plotMutWaves(pdfPath=file.path(plotDir,"wavesall.pdf"),
             variants=list(drivers=0:18, superdrivers=0:7,
                           mutators=NA, passengers=NA))

plotMutWaves(pdfPath=file.path(plotDir,"wavesall_log.pdf"),
             variants=list(drivers=0:18, superdrivers=0:7,
                           mutators=NA, passengers=NA),
             plotY=log10)

## driver passengers ###

plotMutWaves(pdfPath=file.path(plotDir,"wavesall_krkp.pdf"),
             variants=list(drivers=0:10, superdrivers=NA,
                           mutators=NA, passengers=0:20),
             fixedParams=c(driver=0.01, factor=1, 
                           ratenormal=0.00000001, # 1e-8
                           mutator=1000))

plotMutWaves(pdfPath=file.path(plotDir,"wavesall_log_krkp.pdf"),
             variants=list(drivers=0:10, superdrivers=NA,
                           mutators=NA, passengers=0:20),
             fixedParams=c(driver=0.01, factor=1, 
                           ratenormal=0.00000001, # 1e-8
                           mutator=1000),
             plotY=log10)

#############################
sourcef("plotPlentyWaves.R")
#############################

# low superdriver advantage
plotMutWaves(pdfPath=file.path(plotDir,"waves_c1.5.pdf"), 
#              variants=list(#drivers=0:7, superdrivers=0:3,
#                            mutators=NA, passengers=NA),
             fixedParams=c(driver=0.01, factor=1.5, 
                           ratenormal=0.00000001, # 1e-8
                           mutator=1000))

plotMutWaves(pdfPath=file.path(plotDir,"waves_c1.5_log.pdf"), 
             #              variants=list(#drivers=0:7, superdrivers=0:3,
             #                            mutators=NA, passengers=NA),
             fixedParams=c(driver=0.01, factor=1.5, 
                           ratenormal=0.00000001, # 1e-8
                           mutator=1000),
             plotY=log10)

# high superdriver advantage, low driver
plotMutWaves(pdfPath=file.path(plotDir,"waves_c3d0.001.pdf"), 
#              variants=list(#drivers=0:7, superdrivers=0:3,
#                            mutators=NA, passengers=NA),
             fixedParams=c(driver=0.001, factor=3, 
                           ratenormal=0.00000001, # 1e-8
                           mutator=1000))

# high superdriver advantage, high driver
plotMutWaves(pdfPath=file.path(plotDir,"waves_c3d0.1.pdf"), 
#              variants=list(#drivers=0:7, superdrivers=0:3,
#                            mutators=NA, passengers=NA),
             fixedParams=c(driver=0.1, factor=3, 
                           ratenormal=0.00000001, # 1e-8
                           mutator=1000))

# no superdriver advantage, low driver
plotMutWaves(pdfPath=file.path(plotDir,"waves_c1d0.001.pdf"), 
#              variants=list(#drivers=0:7, superdrivers=0:3,
#                            mutators=NA, passengers=NA),
             fixedParams=c(driver=0.001, factor=1, 
                           ratenormal=0.00000001, # 1e-8
                           mutator=1000))

# no superdriver advantage, low driver
plotMutWaves(pdfPath=file.path(plotDir,"waves_c1d0.1.pdf"), 
#              variants=list(#drivers=0:7, superdrivers=0:3,
#                            mutators=NA, passengers=NA),
             fixedParams=c(driver=0.1, factor=1, 
                           ratenormal=0.00000001, # 1e-8
                           mutator=1000))

# low superdriver advantage, low driver
plotMutWaves(pdfPath=file.path(plotDir,"waves_c1d0.001.pdf"), 
#              variants=list(#drivers=0:7, superdrivers=0:3,
#                            mutators=NA, passengers=NA),
             fixedParams=c(driver=0.001, factor=1.5, 
                           ratenormal=0.00000001, # 1e-8
                           mutator=1000))

# low superdriver advantage, high driver
# plotMutWaves(pdfPath=file.path(plotDir,"waves_c1d0.1.pdf"), 
# #              variants=list(#drivers=0:7, superdrivers=0:3,
# #                            mutators=NA, passengers=NA),
#              fixedParams=c(driver=0.1, factor=1.5, 
#                            ratenormal=0.00000001, # 1e-8
#                            mutator=1000))

## example plots
## u=1e-8
plotMutWaves(pdfPath=file.path(plotDir,"waves1a.pdf"))
plotMutWaves(pdfPath=file.path(plotDir,"waves1b.pdf"), isCol="d1", isPch="d2")
plotMutWaves(pdfPath=file.path(plotDir,"waves1c.pdf"), plotY=log10)

plotMutWaves(pdfPath=file.path(plotDir,"waves2a.pdf"), variants=list(drivers=0:7, superdrivers=0:3,
                                                                   mutators=NA, passengers=NA))
plotMutWaves(pdfPath=file.path(plotDir,"waves2b.pdf"), variants=list(drivers=0:7, superdrivers=0:3,
                                                                     mutators=NA, passengers=NA), 
             isCol="d1", isPch="d2")
## u=1e-7
plotMutWaves(pdfPath=file.path(plotDir,"waves3a.pdf"), variants=list(drivers=0:7, superdrivers=0:5,
                                                                     mutators=NA, passengers=NA),
             fixedParams=c(driver=0.01, factor=2, 
                           ratenormal=0.0000001, # 1e-7
                           mutator=1000))
plotMutWaves(pdfPath=file.path(plotDir,"waves3b.pdf"), variants=list(drivers=0:7, superdrivers=0:5,
                                                                     mutators=NA, passengers=NA), 
             isCol="d1", isPch="d2",fixedParams=c(driver=0.01, factor=2, 
                                                  ratenormal=0.0000001, # 1e-7
                                                  mutator=1000))
plotMutWaves(pdfPath=file.path(plotDir,"waves3c.pdf"), variants=list(drivers=0:7, superdrivers=0:5,
                                                                     mutators=NA, passengers=NA),
             fixedParams=c(driver=0.01, factor=2, 
                           ratenormal=0.0000001, # 1e-7
                           mutator=1000), plotY=log10)

##############
### Diversity
##############

#---------------------
### final populations 
#---------------------

### mutation classes ###

plot4mutClassResults.final <- function(case, onPlot, onX, fix, epsil=0.01, parallel=goParallel, 
                                       cores=nCores, xloga=NULL, yloga=NULL, ...) {
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
  
#   browser()
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
                                       epsil=epsil, plotSd=T, xloga=xloga, yloga=yloga, ...)
  #   dev.off()
  objects_numberClasses
}


# finalPops_A_mutator1000driver0.001factor1ratenormal0.00000001_4mutClassesResFinalPops
# finalPops_A_mutator1000driver0.001factor1ratenormal0.00000001_4mutClassesResFinalPops
# goParallel=FALSE
# plot4mutClassResults.final("A_", fix=c("ratenormal"=0.00000001), onPlot=c("factor","mutator"), 
#                            onX="driver", epsil=0.001)
# plot4mutClassResults.final("B_", fix=c("ratenormal"=0.00000001), onPlot=c("factor","mutator"), 
#                            onX="driver", epsil=0.001)
# plot4mutClassResults.final("B_", fix=c("ratenormal"=0.00000001), onPlot=c("factor","mutator"), 
#                            onX="driver", epsil=0.001, xloga=log)
# plot4mutClassResults.final("C_", fix=c("ratenormal"=0.00000001), onPlot=c("mutator","factor"), 
#                            onX="driver", epsil=0.001)
# plot4mutClassResults.final("D_", fix=c("ratenormal"=0.00000001), onPlot=c("mutator","factor"), 
#                            onX="driver", epsil=0.001)

plotn4mutClasses42plot <- function() {
  cases <- sapply(LETTERS[1:4], function(x) paste0(x,"_"))
  casesPlotnames <- sapply(cases, function(x) sprintf("%s4mutClassesResFinalPops",x))
  plotFun4thisParamValues(cases=cases, casesPlotnames=casesPlotnames, plot4mutClassResults.final)
}

######################################
### plot the four mutation classes ###
plotn4mutClasses42plot()
######################################

### average number of (super)drivers ###

plotOneClones_highestFreq.final <- function(case, onPlot, onX, fix, epsil=0.01, parallel=goParallel, 
                                            cores=nCores, xloga=NULL, yloga=NULL, ...) {  
  objects <- tmp_getObjects(sprintf("%s.*oneClonesFinalPopsRes.*\\.rda",case), parallel=parallel, cores=cores)
  
  #   print("the oneConesFinalRes files")
  #   print(fileNames)
  
  if (length(objects) == 0) {
    print(sprintf("length objects 0, with oneClones of case %s. nothing to plot!",case))
#     browser()
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
                                             xloga=xloga, yloga=yloga, 
                                             legendApprox="bottomright", ...)
  plot_varySome.oneClonesMutations.finalPops(objects_highest_drivers, 
                                             onPlot=onPlot, onX=onX, fix=fix, 
                                             whichMutation="drivers", 
                                             epsil=epsil, plotSd=T,
                                             xloga=xloga, yloga=yloga,
                                             legendApprox="bottomright", ...)
  plot_varySome.oneClonesMutations.finalPops(objects_highest_mutators, 
                                             onPlot=onPlot, onX=onX, fix=fix, 
                                             whichMutation="mutators", 
                                             epsil=epsil, plotSd=FALSE,
                                             xloga=xloga, yloga=yloga,
                                             legendApprox="bottomright", ...)
  #   dev.off()
  objects_highest
}

# tmp11 <- plotOneClones_highestFreq.final("A_", fix=c("ratenormal"=0.00000001), onPlot=c("factor","mutator"), 
#                                          onX="driver", epsil=0.001)
# tmp11 <- plotOneClones_highestFreq.final("A_", fix=c("ratenormal"=0.00000001), onPlot=c("factor","mutator"), 
#                                          onX="driver", epsil=0.001)
# tmp22 <- plotOneClones_highestFreq.final("B_", fix=c("ratenormal"=0.00000001), onPlot=c("factor","mutator"), 
#                                          onX="driver", epsil=0.001)
# tmp33 <- plotOneClones_highestFreq.final("C_", fix=c("ratenormal"=0.00000001), onPlot=c("mutator","factor"), 
#                                          onX="driver", epsil=0.001)
# tmp44 <- plotOneClones_highestFreq.final("D_", fix=c("ratenormal"=0.00000001), onPlot=c("mutator","factor"), 
#                                          onX="driver", epsil=0.001)

plotOneCloness42plot <- function() {
  cases <- sapply(LETTERS[1:4], function(x) paste0(x,"_"))
  casesPlotnames <- sapply(cases, function(x) sprintf("%soneClonesFinalPopsRes_leading",x))
  plotFun4thisParamValues(cases=cases, casesPlotnames=casesPlotnames, plotOneClones_highestFreq.final)
}

###################################
### plot the leading one clones ###
plotOneCloness42plot()
###################################


########################
#### leading clone #####

plotOneClones.bar <- function(cases, N=50) {
  require(ggplot2)
  
  #' QUICK AND DIRTY!!!
  #' @param cases vector of names, such as "B_mutator1000driver0.01factor1.5ratenormal0.00000001"
     #' Note please that the case should be included.
     getLeadingCloneFrequencies <- function(cases) {
       
       res <- sapply(cases, function(x) {
         obj <- tryCatch(
           loadIfNot(file.path(intermediateDir,sprintf("finalPops_%s_oneClonesFinalPops_highest.rda",x))),
                  error=function(e) return(NA))
         if (!is.na(obj)) {
           vals <- sapply(obj, function(x) x$mean.count)
           mean <- mean(vals)
           sd <- sd(vals)
           se <- sd/sqrt(N)
           c(mean=mean,sd=sd,se=se)
         } else return(c(mean=NA,sd=NA,se=NA))        
       })
       
#        res <- Filter(function(x) !is.na(x), res)
       
       colnames(res) <- sub("[ABCD]_(.+)","\\1",getShortStrings(colnames(res)))
       colnames(res) <- sub("0.00000001", "1e-8", colnames(res))
       colnames(res) <- sub("0.0000001", "1e-7", colnames(res))  
       colnres <- colnames(res)
#        colnames(res) <- sapply(1:length(colnres), function(i) sprintf("%s) %s", letters[i], colnres[i]))
       colnames(res) <- sapply(1:length(colnres), function(i) sprintf("%s) %s", letters[i], colnres[i]))
        res <- res[,rev(colnames(res))]
	t(res)
     }
     
     res <- getLeadingCloneFrequencies(cases)
#      res <- cbind(comb=rownames(res),res)
     
#   xx <- factor(rownames(res))
#   browser()
#       res[is.na(res)] <- 0
#   rownres <- rownames(res)
#   rownames(res) <- sapply(1:length(rownres), function(i) sprintf("%s) %s", letters[i], rownres[i]))
  print(res)
     ggplot(as.data.frame(res), aes(x =factor(rownames(res)), y = mean), environment=environment() ) +
       geom_bar(position = position_dodge(), stat="identity") +
       geom_errorbar(aes(ymin=mean-se, ymax=mean+se)) +
       xlab("Parameter combination") + ylab("Mean frequency [%]") +
       theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
       coord_flip()
}

oneClonesForFreq <- c(
  "A_mutator1000driver0.001factor1ratenormal0.00000001",
  "A_mutator1000driver0.001factor1ratenormal0.0000001",
  "A_mutator1000driver0.01factor1ratenormal0.00000001",
  "A_mutator1000driver0.01factor1ratenormal0.0000001",
  "A_mutator1000driver0.1factor1ratenormal0.00000001",
  "A_mutator1000driver0.1factor1ratenormal0.0000001",
  #
  "B_mutator1000driver0.001factor2ratenormal0.00000001",
  "B_mutator1000driver0.01factor2ratenormal0.00000001",
  "B_mutator1000driver0.1factor2ratenormal0.00000001",
  #
  "C_mutator1driver0.001factor2ratenormal0.00000001",
  "C_mutator1driver0.01factor2ratenormal0.00000001",
  "C_mutator1driver0.1factor2ratenormal0.00000001"
  )

freqBars <- plotOneClones.bar(oneClonesForFreq)
ggsave(file.path(plotDir,"barplot_oneClonesFreq.pdf"), plot=freqBars)

############################
### most frequent clones ###

tmp_getPops <- function(fullNameCases) {
  pops <- lapply(fullNameCases, function(x) {
    fileName <- file.path(rdataDir,paste0("finalPops_",x,".rda"))
    print(paste0("load file ", fileName))
    loadIfNot(fileName)
  })
#   names(pops) <- sub("0.00000001", "1e-08",fullNameCases)
  names(pops) <- c("Model A", "Model B", "Model C", "Model D")
  pops
}

somepops <- tmp_getPops(c("A_mutator1000driver0.01factor1ratenormal0.00000001",
                      "B_mutator1000driver0.001factor2ratenormal0.00000001",
                      "C_mutator5driver0.01factor2ratenormal0.00000001",
                      "D_mutator5driver0.01factor1ratenormal0.00000001"))

somepops_up <- getParts2(somepops, "upper", 50)
somepops_up_means <- getMeansOf2(somepops_up)
tmp_somecompareMeans <- list(A=somepops_up_means)
pdf(file.path(plotDir,"someMostFrequent.pdf"))
scatterPlot(tmp_somecompareMeans, columns=c("count"), histType="highLow", plotNamesAs=TRUE)
dev.off()

