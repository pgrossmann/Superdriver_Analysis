##################################################################################
###
### Here are the wave plots for publication
###
### Patrick Grossmann, December 13, 2016
###
##################################################################################

#######################################################################################
####### --------------------------- Parameter set 1 --------------------------- #######
fixed <- c(driver=0.01, ratenormal=0.00000001, mutator=1000)
variants <- list(drivers=0:10, superdrivers=0:5, mutators=NA, passengers=NA)

# plotWaveStats()

temp_cases <- c("B_mutator1000driver0.01factor1.1ratenormal0.00000001_4mutClasses",
                "B_mutator1000driver0.01factor1.3ratenormal0.00000001_4mutClasses",
                "B_mutator1000driver0.01factor1.5ratenormal0.00000001_4mutClasses",
                "B_mutator1000driver0.01factor2.6ratenormal0.00000001_4mutClasses",
                "B_mutator1000driver0.01factor2.8ratenormal0.00000001_4mutClasses",
                "B_mutator1000driver0.01factor3ratenormal0.00000001_4mutClasses")

plotQDistributions(temp_cases, plotAsterisk=F)

getQuadraticDistribution(case="B_mutator1000driver0.01factor1.1ratenormal0.00000001_4mutClasses")

#--------------------------------------------------------------------#
### Driver 0.01, superdriver factor 2 (Figure Waves per Replicate) ###
### This is Edit August 18, 2017 

repNums <- 50:1

## All waves but only one replicate
for (repNum in repNums) {
  
  # plotMutWaves(pdfPath=file.path(plotDir, sprintf("waves_c2_rep%s_log.pdf", repNum)), 
  #              fixedParams=c(fixed, factor=1.1), variants=variants, plotY=log10,
  #              whichReplicate=repNum)
  # plotMutWaves(pdfPath=file.path(plotDir, sprintf("waves_c2_rep%s_log.pdf", repNum)),
  #              fixedParams=c(fixed, factor=1.3), variants=variants, plotY=log10,
  #              whichReplicate=repNum)
  # plotMutWaves(pdfPath=file.path(plotDir, sprintf("waves_c2_rep%s_log.pdf", repNum)),
  #              fixedParams=c(fixed, factor=1.5), variants=variants, plotY=log10,
  #              whichReplicate=repNum)
  # plotMutWaves(pdfPath=file.path(plotDir, sprintf("waves_c2_rep%s_log.pdf", repNum)),
  #              fixedParams=c(fixed, factor=2.6), variants=variants, plotY=log10,
  #              whichReplicate=repNum)
  # plotMutWaves(pdfPath=file.path(plotDir, sprintf("waves_c2_rep%s_log.pdf", repNum)),
  #              fixedParams=c(fixed, factor=2.8), variants=variants, plotY=log10,
  #              whichReplicate=repNum)
  # plotMutWaves(pdfPath=file.path(plotDir, sprintf("waves_c2_rep%s_log.pdf", repNum)),
  #              fixedParams=c(fixed, factor=3), variants=variants, plotY=log10,
  #              whichReplicate=repNum)
}

#--------------------------------------#
### Vary superdriver (Figure Shift) ###

# Low superdriver advantages
cc <- c(1.1, 1.3, 1.5)
sapply(cc, function(c) plotMutWaves(pdfPath=file.path(plotDir,sprintf("waves_c%s_log.pdf", c)), 
                                    fixedParams=c(fixed, factor=c), 
                                    variants=variants, plotY=log10))

# High superdriver advantages
cc <- c(2.6, 2.8, 3)
sapply(cc, function(c) plotMutWaves(pdfPath=file.path(plotDir,sprintf("waves_c%s_log.pdf", c)),
                                    fixedParams=c(fixed, factor=c), 
                                    variants=variants, plotY=log10))

#---------------------------------------------------#
### Vary superdriver (Figure Shift) for many reps ###
for (repNum in repNums) {
  
  # Low superdriver advantages
  cc <- c(1.1, 1.3, 1.5)
  sapply(cc, function(c) plotMutWaves(pdfPath=file.path(plotDir,sprintf("waves_c%s_rep%s_log.pdf", c, repNum)), 
                                      fixedParams=c(fixed, factor=c), 
                                      variants=variants, plotY=log10,
                                      whichReplicate=repNum))
  
  # High superdriver advantages
  cc <- c(2.6, 2.8, 3)
  sapply(cc, function(c) plotMutWaves(pdfPath=file.path(plotDir,sprintf("waves_c%s_rep%s_log.pdf", c, repNum)),
                                      fixedParams=c(fixed, factor=c), 
                                      variants=variants, plotY=log10,
                                      whichReplicate=repNum))
}

#--------------------------------------#
### Vary superdriver (Figure Shift) ###

# Low superdriver advantages
cc <- c(1.1, 1.3, 1.5)
sapply(cc, function(c) plotMutWaves(pdfPath=file.path(plotDir,sprintf("waves_c%s_log.pdf", c)), 
                                    fixedParams=c(fixed, factor=c), 
                                    variants=variants, plotY=log10))

# High superdriver advantages
cc <- c(2.6, 2.8, 3)
sapply(cc, function(c) plotMutWaves(pdfPath=file.path(plotDir,sprintf("waves_c%s_log.pdf", c)),
                                    fixedParams=c(fixed, factor=c), 
                                    variants=variants, plotY=log10))


#---------------------------------------------------------#
### Driver 0.01, superdriver factor 2 (Figure Overview) ###

## All waves
plotMutWaves(pdfPath=file.path(plotDir, "waves_c2.pdf"), 
             fixedParams=c(fixed, factor=2), variants=variants)
plotMutWaves(pdfPath=file.path(plotDir, "waves_c2_log.pdf"), 
             fixedParams=c(fixed, factor=2), variants=variants, plotY=log10)

## All waves, supplements
dd <- c(0.005, 0.02, 0.03, 0.04, 0.05)
sapply(dd, function(dr) plotMutWaves(pdfPath=file.path(plotDir, sprintf("waves_dr%s_c2.pdf", dr)),
                                     fixedParams=c(driver=dr, ratenormal=0.00000001, mutator=1000, factor=2),
                                     variants=variants, plotY=log10))

## Filter waves
# for filtering 
tmpListS0 <- lapply(0:100, function(x) c(x,0)) # superdriver0
tmpListS1 <- lapply(0:100, function(x) c(x,1)) # superdriver1
tmpListS2 <- lapply(0:100, function(x) c(x,2)) # superdriver2
tmpListS3 <- lapply(0:100, function(x) c(x,3)) # superdriver3
tmpListS4 <- lapply(0:100, function(x) c(x,4)) # superdriver4
tmpListS5 <- lapply(0:100, function(x) c(x,5)) # superdriver5

tmpListD0 <- lapply(0:10, function(x) c(0,x)) # driver0
tmpListD1 <- lapply(0:10, function(x) c(1,x)) # driver1
tmpListD2 <- lapply(0:10, function(x) c(2,x)) # driver2
tmpListD3 <- lapply(0:10, function(x) c(3,x)) # driver3
tmpListD4 <- lapply(0:10, function(x) c(4,x)) # driver4
tmpListD5 <- lapply(0:10, function(x) c(5,x)) # driver5
tmpListD6 <- lapply(0:10, function(x) c(6,x)) # driver6

# filtered conditional on superdrivers
plotMutWaves(pdfPath=file.path(plotDir,"waves_c2_filtered_superdriver0_log.pdf"), 
             #              variants=list(#drivers=0:7, superdrivers=0:3,
             #                            mutators=NA, passengers=NA),
             fixedParams=c(driver=0.01, factor=2, 
                           ratenormal=0.00000001, # 1e-8
                           mutator=1000),
             filter2dVec_in=tmpListS0, plotY=log10)
plotMutWaves(pdfPath=file.path(plotDir,"waves_c2_filtered_superdriver2_log.pdf"), 
             #              variants=list(#drivers=0:7, superdrivers=0:3,
             #                            mutators=NA, passengers=NA),
             fixedParams=c(driver=0.01, factor=2, 
                           ratenormal=0.00000001, # 1e-8
                           mutator=1000),
             filter2dVec_in=tmpListS2, plotY=log10)
plotMutWaves(pdfPath=file.path(plotDir,"waves_c2_filtered_superdriver4_log.pdf"), 
             #              variants=list(#drivers=0:7, superdrivers=0:3,
             #                            mutators=NA, passengers=NA),
             fixedParams=c(driver=0.01, factor=2, 
                           ratenormal=0.00000001, # 1e-8
                           mutator=1000),
             filter2dVec_in=tmpListS4, plotY=log10)

# filtered conditional on drivers
plotMutWaves(pdfPath=file.path(plotDir,"waves_c2_filtered_driver0_log.pdf"), 
             #              variants=list(#drivers=0:7, superdrivers=0:3,
             #                            mutators=NA, passengers=NA),
             fixedParams=c(driver=0.01, factor=2, 
                           ratenormal=0.00000001, # 1e-8
                           mutator=1000),
             filter2dVec_in=tmpListD0, plotY=log10)
plotMutWaves(pdfPath=file.path(plotDir,"waves_c2_filtered_driver2_log.pdf"), 
             #              variants=list(#drivers=0:7, superdrivers=0:3,
             #                            mutators=NA, passengers=NA),
             fixedParams=c(driver=0.01, factor=2, 
                           ratenormal=0.00000001, # 1e-8
                           mutator=1000),
             filter2dVec_in=tmpListD2, plotY=log10)
plotMutWaves(pdfPath=file.path(plotDir,"waves_c2_filtered_driver4_log.pdf"), 
             #              variants=list(#drivers=0:7, superdrivers=0:3,
             #                            mutators=NA, passengers=NA),
             fixedParams=c(driver=0.01, factor=2, 
                           ratenormal=0.00000001, # 1e-8
                           mutator=1000),
             filter2dVec_in=tmpListD4, plotY=log10)


#----------------------------------------------------------#
### Driver 0.01, superdriver factor 2 (Figure Diversity) ###
source(file.path(scriptDir, "analyseMore_functions.R"))

for (k in 1:3) {
  for (l in 1:3) {
    pdf(file.path(plotDir, sprintf("B_waitingTimeResk%sl%s.pdf",k,l)))
    # plotWaitingTimesResults_2d("B_.+factor(1.5|2|2.5|3)r", fix=c("mutator"=1000,"ratenormal"=0.00000001),
    #                            onPlot=c("factor"), onX="driver", epsil=0.001, kCancer=c(kr=k,ks=l))
    dev.off()
  }
}

#----------------------------------------------------------#
### Driver 0.01, superdriver factor 2 (Figure Diversity) ###
source(file.path(scriptDir, "analyseMore_functions.R"))

#plotDiversity()