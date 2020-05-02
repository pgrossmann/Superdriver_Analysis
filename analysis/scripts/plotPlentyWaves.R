## for filtering ##
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
###

# 2-fold superdriver advantage
plotMutWaves(pdfPath=file.path(plotDir,"waves_c2.pdf"), 
             #              variants=list(#drivers=0:7, superdrivers=0:3,
             #                            mutators=NA, passengers=NA),
             fixedParams=c(driver=0.01, factor=2, 
                           ratenormal=0.00000001, # 1e-8
                           mutator=1000))
plotMutWaves(pdfPath=file.path(plotDir,"waves_c2_log.pdf"), 
             #              variants=list(#drivers=0:7, superdrivers=0:3,
             #                            mutators=NA, passengers=NA),
             fixedParams=c(driver=0.01, factor=2, 
                           ratenormal=0.00000001, # 1e-8
                           mutator=1000),
             plotY=log10)
## filtered
plotMutWaves(pdfPath=file.path(plotDir,"waves_c2_filtered_superdriver0.pdf"), 
             #              variants=list(#drivers=0:7, superdrivers=0:3,
             #                            mutators=NA, passengers=NA),
             fixedParams=c(driver=0.01, factor=2, 
                           ratenormal=0.00000001, # 1e-8
                           mutator=1000),
             filter2dVec_in=tmpListS0)
plotMutWaves(pdfPath=file.path(plotDir,"waves_c2_filtered_superdriver1.pdf"), 
             #              variants=list(#drivers=0:7, superdrivers=0:3,
             #                            mutators=NA, passengers=NA),
             fixedParams=c(driver=0.01, factor=2, 
                           ratenormal=0.00000001, # 1e-8
                           mutator=1000),
             filter2dVec_in=tmpListS1)
plotMutWaves(pdfPath=file.path(plotDir,"waves_c2_filtered_superdriver4.pdf"), 
             #              variants=list(#drivers=0:7, superdrivers=0:3,
             #                            mutators=NA, passengers=NA),
             fixedParams=c(driver=0.01, factor=2, 
                           ratenormal=0.00000001, # 1e-8
                           mutator=1000),
             filter2dVec_in=tmpListS4)
plotMutWaves(pdfPath=file.path(plotDir,"waves_c2_filtered_superdriver2.pdf"), 
             #              variants=list(#drivers=0:7, superdrivers=0:3,
             #                            mutators=NA, passengers=NA),
             fixedParams=c(driver=0.01, factor=2, 
                           ratenormal=0.00000001, # 1e-8
                           mutator=1000),
             filter2dVec_in=tmpListS2)
plotMutWaves(pdfPath=file.path(plotDir,"waves_c2_filtered_superdriver3.pdf"), 
             #              variants=list(#drivers=0:7, superdrivers=0:3,
             #                            mutators=NA, passengers=NA),
             fixedParams=c(driver=0.01, factor=2, 
                           ratenormal=0.00000001, # 1e-8
                           mutator=1000),
             filter2dVec_in=tmpListS3)
plotMutWaves(pdfPath=file.path(plotDir,"waves_c2_filtered_superdriver5.pdf"), 
             #              variants=list(#drivers=0:7, superdrivers=0:3,
             #                            mutators=NA, passengers=NA),
             fixedParams=c(driver=0.01, factor=2, 
                           ratenormal=0.00000001, # 1e-8
                           mutator=1000),
             filter2dVec_in=tmpListS5)
#
plotMutWaves(pdfPath=file.path(plotDir,"waves_c2_filtered_superdriver0_log.pdf"), 
             #              variants=list(#drivers=0:7, superdrivers=0:3,
             #                            mutators=NA, passengers=NA),
             fixedParams=c(driver=0.01, factor=2, 
                           ratenormal=0.00000001, # 1e-8
                           mutator=1000),
             filter2dVec_in=tmpListS0,
             plotY=log10)
plotMutWaves(pdfPath=file.path(plotDir,"waves_c2_filtered_superdriver1_log.pdf"), 
             #              variants=list(#drivers=0:7, superdrivers=0:3,
             #                            mutators=NA, passengers=NA),
             fixedParams=c(driver=0.01, factor=2, 
                           ratenormal=0.00000001, # 1e-8
                           mutator=1000),
             filter2dVec_in=tmpListS1,
             plotY=log10)
plotMutWaves(pdfPath=file.path(plotDir,"waves_c2_filtered_superdriver4_log.pdf"), 
             #              variants=list(#drivers=0:7, superdrivers=0:3,
             #                            mutators=NA, passengers=NA),
             fixedParams=c(driver=0.01, factor=2, 
                           ratenormal=0.00000001, # 1e-8
                           mutator=1000),
             filter2dVec_in=tmpListS4,
             plotY=log10)
plotMutWaves(pdfPath=file.path(plotDir,"waves_c2_filtered_superdriver2_log.pdf"), 
             #              variants=list(#drivers=0:7, superdrivers=0:3,
             #                            mutators=NA, passengers=NA),
             fixedParams=c(driver=0.01, factor=2, 
                           ratenormal=0.00000001, # 1e-8
                           mutator=1000),
             filter2dVec_in=tmpListS2,
             plotY=log10)
plotMutWaves(pdfPath=file.path(plotDir,"waves_c2_filtered_superdriver3_log.pdf"), 
             #              variants=list(#drivers=0:7, superdrivers=0:3,
             #                            mutators=NA, passengers=NA),
             fixedParams=c(driver=0.01, factor=2, 
                           ratenormal=0.00000001, # 1e-8
                           mutator=1000),
             filter2dVec_in=tmpListS3,
             plotY=log10)
#
plotMutWaves(pdfPath=file.path(plotDir,"waves_c2_filtered_superdriver3superdriver4.pdf"), 
             #              variants=list(#drivers=0:7, superdrivers=0:3,
             #                            mutators=NA, passengers=NA),
             fixedParams=c(driver=0.01, factor=2, 
                           ratenormal=0.00000001, # 1e-8
                           mutator=1000),
             filter2dVec_in=c(tmpListS3,tmpListS4))
plotMutWaves(pdfPath=file.path(plotDir,"waves_c2_filtered_superdriver3superdriver4_log.pdf"), 
             #              variants=list(#drivers=0:7, superdrivers=0:3,
             #                            mutators=NA, passengers=NA),
             fixedParams=c(driver=0.01, factor=2, 
                           ratenormal=0.00000001, # 1e-8
                           mutator=1000),
             filter2dVec_in=c(tmpListS3,tmpListS4),
             plotY=log10)
plotMutWaves(pdfPath=file.path(plotDir,"waves_c2_filtered_superdriver1superdriver4_log.pdf"), 
             #              variants=list(#drivers=0:7, superdrivers=0:3,
             #                            mutators=NA, passengers=NA),
             fixedParams=c(driver=0.01, factor=2, 
                           ratenormal=0.00000001, # 1e-8
                           mutator=1000),
             filter2dVec_in=c(tmpListS1,tmpListS4),
             plotY=log10)

#
plotMutWaves(pdfPath=file.path(plotDir,"waves_c2_filtered_driver0.pdf"), 
             #              variants=list(#drivers=0:7, superdrivers=0:3,
             #                            mutators=NA, passengers=NA),
             fixedParams=c(driver=0.01, factor=2, 
                           ratenormal=0.00000001, # 1e-8
                           mutator=1000),
             filter2dVec_in=tmpListD0)
plotMutWaves(pdfPath=file.path(plotDir,"waves_c2_filtered_driver1.pdf"), 
             #              variants=list(#drivers=0:7, superdrivers=0:3,
             #                            mutators=NA, passengers=NA),
             fixedParams=c(driver=0.01, factor=2, 
                           ratenormal=0.00000001, # 1e-8
                           mutator=1000),
             filter2dVec_in=tmpListD1)
plotMutWaves(pdfPath=file.path(plotDir,"waves_c2_filtered_driver2.pdf"), 
             #              variants=list(#drivers=0:7, superdrivers=0:3,
             #                            mutators=NA, passengers=NA),
             fixedParams=c(driver=0.01, factor=2, 
                           ratenormal=0.00000001, # 1e-8
                           mutator=1000),
             filter2dVec_in=tmpListD2)
plotMutWaves(pdfPath=file.path(plotDir,"waves_c2_filtered_driver3.pdf"), 
             #              variants=list(#drivers=0:7, superdrivers=0:3,
             #                            mutators=NA, passengers=NA),
             fixedParams=c(driver=0.01, factor=2, 
                           ratenormal=0.00000001, # 1e-8
                           mutator=1000),
             filter2dVec_in=tmpListD3)
plotMutWaves(pdfPath=file.path(plotDir,"waves_c2_filtered_driver4.pdf"), 
             #              variants=list(#drivers=0:7, superdrivers=0:3,
             #                            mutators=NA, passengers=NA),
             fixedParams=c(driver=0.01, factor=2, 
                           ratenormal=0.00000001, # 1e-8
                           mutator=1000),
             filter2dVec_in=tmpListD4)
plotMutWaves(pdfPath=file.path(plotDir,"waves_c2_filtered_driver5.pdf"), 
             #              variants=list(#drivers=0:7, superdrivers=0:3,
             #                            mutators=NA, passengers=NA),
             fixedParams=c(driver=0.01, factor=2, 
                           ratenormal=0.00000001, # 1e-8
                           mutator=1000),
             filter2dVec_in=tmpListD5)
#
plotMutWaves(pdfPath=file.path(plotDir,"waves_c2_filtered_driver0_log.pdf"), 
             #              variants=list(#drivers=0:7, superdrivers=0:3,
             #                            mutators=NA, passengers=NA),
             fixedParams=c(driver=0.01, factor=2, 
                           ratenormal=0.00000001, # 1e-8
                           mutator=1000),
             filter2dVec_in=tmpListD0,
             plotY=log10)
plotMutWaves(pdfPath=file.path(plotDir,"waves_c2_filtered_driver1_log.pdf"), 
             #              variants=list(#drivers=0:7, superdrivers=0:3,
             #                            mutators=NA, passengers=NA),
             fixedParams=c(driver=0.01, factor=2, 
                           ratenormal=0.00000001, # 1e-8
                           mutator=1000),
             filter2dVec_in=tmpListD1,
             plotY=log10)
plotMutWaves(pdfPath=file.path(plotDir,"waves_c2_filtered_driver2_log.pdf"), 
             #              variants=list(#drivers=0:7, superdrivers=0:3,
             #                            mutators=NA, passengers=NA),
             fixedParams=c(driver=0.01, factor=2, 
                           ratenormal=0.00000001, # 1e-8
                           mutator=1000),
             filter2dVec_in=tmpListD2,
             plotY=log10)
plotMutWaves(pdfPath=file.path(plotDir,"waves_c2_filtered_driver3_log.pdf"), 
             #              variants=list(#drivers=0:7, superdrivers=0:3,
             #                            mutators=NA, passengers=NA),
             fixedParams=c(driver=0.01, factor=2, 
                           ratenormal=0.00000001, # 1e-8
                           mutator=1000),
             filter2dVec_in=tmpListD3,
             plotY=log10)
plotMutWaves(pdfPath=file.path(plotDir,"waves_c2_filtered_driver4_log.pdf"), 
             #              variants=list(#drivers=0:7, superdrivers=0:3,
             #                            mutators=NA, passengers=NA),
             fixedParams=c(driver=0.01, factor=2, 
                           ratenormal=0.00000001, # 1e-8
                           mutator=1000),
             filter2dVec_in=tmpListD4,
             plotY=log10)
#
plotMutWaves(pdfPath=file.path(plotDir,"waves_c2_filtered_driver1_superdriver4.pdf"), 
             #              variants=list(#drivers=0:7, superdrivers=0:3,
             #                            mutators=NA, passengers=NA),
             fixedParams=c(driver=0.01, factor=2, 
                           ratenormal=0.00000001, # 1e-8
                           mutator=1000),
             filter2dVec_in=c(tmpListS4,tmpListD1))
# plot them log
s0 <-
  plotMutWaves(pdfPath=file.path(plotDir,"waves_c2_filtered_superdriver0_log.pdf"), 
               #              variants=list(#drivers=0:7, superdrivers=0:3,
               #                            mutators=NA, passengers=NA),
               fixedParams=c(driver=0.01, factor=2, 
                             ratenormal=0.00000001, # 1e-8
                             mutator=1000),
               filter2dVec_in=tmpListS0,
               plotY=log10)
s1 <-
  plotMutWaves(pdfPath=file.path(plotDir,"waves_c2_filtered_superdriver1_log.pdf"), 
               #              variants=list(#drivers=0:7, superdrivers=0:3,
               #                            mutators=NA, passengers=NA),
               fixedParams=c(driver=0.01, factor=2, 
                             ratenormal=0.00000001, # 1e-8
                             mutator=1000),
               filter2dVec_in=tmpListS1,
               plotY=log10)

s2 <-
  plotMutWaves(pdfPath=file.path(plotDir,"waves_c2_filtered_superdriver2_log.pdf"), 
               #              variants=list(#drivers=0:7, superdrivers=0:3,
               #                            mutators=NA, passengers=NA),
               fixedParams=c(driver=0.01, factor=2, 
                             ratenormal=0.00000001, # 1e-8
                             mutator=1000),
               filter2dVec_in=tmpListS2,
               plotY=log10)
s3 <-
  plotMutWaves(pdfPath=file.path(plotDir,"waves_c2_filtered_superdriver3_log.pdf"), 
               #              variants=list(#drivers=0:7, superdrivers=0:3,
               #                            mutators=NA, passengers=NA),
               fixedParams=c(driver=0.01, factor=2, 
                             ratenormal=0.00000001, # 1e-8
                             mutator=1000),
               filter2dVec_in=tmpListS3,
               plotY=log10)
s4 <-
  plotMutWaves(pdfPath=file.path(plotDir,"waves_c2_filtered_superdriver4_log.pdf"), 
               #              variants=list(#drivers=0:7, superdrivers=0:3,
               #                            mutators=NA, passengers=NA),
               fixedParams=c(driver=0.01, factor=2, 
                             ratenormal=0.00000001, # 1e-8
                             mutator=1000),
               filter2dVec_in=tmpListS4,
               plotY=log10)
s5 <-
  plotMutWaves(pdfPath=file.path(plotDir,"waves_c2_filtered_superdriver5_log.pdf"), 
               #              variants=list(#drivers=0:7, superdrivers=0:3,
               #                            mutators=NA, passengers=NA),
               fixedParams=c(driver=0.01, factor=2, 
                             ratenormal=0.00000001, # 1e-8
                             mutator=1000),
               filter2dVec_in=tmpListS5,
               plotY=log10)
#
s0_left <-
  plotMutWaves(pdfPath=file.path(plotDir,"waves_c2_filtered_superdriver0_log.pdf"), 
               #              variants=list(#drivers=0:7, superdrivers=0:3,
               #                            mutators=NA, passengers=NA),
               fixedParams=c(driver=0.01, factor=2, 
                             ratenormal=0.00000001, # 1e-8
                             mutator=1000),
               filter2dVec_in=tmpListS0,
               plotY=log10,
               midFunX=min)
s1_left <-
  plotMutWaves(pdfPath=file.path(plotDir,"waves_c2_filtered_superdriver1_log.pdf"), 
               #              variants=list(#drivers=0:7, superdrivers=0:3,
               #                            mutators=NA, passengers=NA),
               fixedParams=c(driver=0.01, factor=2, 
                             ratenormal=0.00000001, # 1e-8
                             mutator=1000),
               filter2dVec_in=tmpListS1,
               plotY=log10,
               midFunX=min)

s2_left <-
  plotMutWaves(pdfPath=file.path(plotDir,"waves_c2_filtered_superdriver2_log.pdf"), 
               #              variants=list(#drivers=0:7, superdrivers=0:3,
               #                            mutators=NA, passengers=NA),
               fixedParams=c(driver=0.01, factor=2, 
                             ratenormal=0.00000001, # 1e-8
                             mutator=1000),
               filter2dVec_in=tmpListS2,
               plotY=log10,
               midFunX=min)
s3_left <-
  plotMutWaves(pdfPath=file.path(plotDir,"waves_c2_filtered_superdriver3_log.pdf"), 
               #              variants=list(#drivers=0:7, superdrivers=0:3,
               #                            mutators=NA, passengers=NA),
               fixedParams=c(driver=0.01, factor=2, 
                             ratenormal=0.00000001, # 1e-8
                             mutator=1000),
               filter2dVec_in=tmpListS3,
               plotY=log10,
               midFunX=min)
s4_left <-
  plotMutWaves(pdfPath=file.path(plotDir,"waves_c2_filtered_superdriver4_log.pdf"), 
               #              variants=list(#drivers=0:7, superdrivers=0:3,
               #                            mutators=NA, passengers=NA),
               fixedParams=c(driver=0.01, factor=2, 
                             ratenormal=0.00000001, # 1e-8
                             mutator=1000),
               filter2dVec_in=tmpListS4,
               plotY=log10,
               midFunX=min)
s5_left <-
  plotMutWaves(pdfPath=file.path(plotDir,"waves_c2_filtered_superdriver5_log.pdf"), 
               #              variants=list(#drivers=0:7, superdrivers=0:3,
               #                            mutators=NA, passengers=NA),
               fixedParams=c(driver=0.01, factor=2, 
                             ratenormal=0.00000001, # 1e-8
                             mutator=1000),
               filter2dVec_in=tmpListS5,
               plotY=log10,
               midFunX=min)
##
plotMutWaves(pdfPath=file.path(plotDir,"waves_c2_filtered_superdriver3superdriver4_log.pdf"), 
             #              variants=list(#drivers=0:7, superdrivers=0:3,
             #                            mutators=NA, passengers=NA),
             fixedParams=c(driver=0.01, factor=2, 
                           ratenormal=0.00000001, # 1e-8
                           mutator=1000),
             filter2dVec_in=c(tmpListS3,tmpListS4),
             plotY=log10)
plotMutWaves(pdfPath=file.path(plotDir,"waves_c2_filtered_superdriver2superdriver3_log.pdf"), 
             #              variants=list(#drivers=0:7, superdrivers=0:3,
             #                            mutators=NA, passengers=NA),
             fixedParams=c(driver=0.01, factor=2, 
                           ratenormal=0.00000001, # 1e-8
                           mutator=1000),
             filter2dVec_in=c(tmpListS2,tmpListS3),
             plotY=log10)
plotMutWaves(pdfPath=file.path(plotDir,"waves_c2_filtered_superdriver4superdriver5_log.pdf"), 
             #              variants=list(#drivers=0:7, superdrivers=0:3,
             #                            mutators=NA, passengers=NA),
             fixedParams=c(driver=0.01, factor=2, 
                           ratenormal=0.00000001, # 1e-8
                           mutator=1000),
             filter2dVec_in=c(tmpListS4,tmpListS5),
             plotY=log10)
#
d0 <-
  plotMutWaves(pdfPath=file.path(plotDir,"waves_c2_filtered_driver0_log.pdf"), 
               variants=list(drivers=0:5, superdrivers=0:7,
                             mutators=NA, passengers=NA),
               fixedParams=c(driver=0.01, factor=2, 
                             ratenormal=0.00000001, # 1e-8
                             mutator=1000),
               filter2dVec_in=tmpListD0,
               plotY=log10,
               plotMidX=FALSE)
d1 <-
  plotMutWaves(pdfPath=file.path(plotDir,"waves_c2_filtered_driver1_log.pdf"), 
               variants=list(drivers=0:5, superdrivers=0:7,
                             mutators=NA, passengers=NA),
               fixedParams=c(driver=0.01, factor=2, 
                             ratenormal=0.00000001, # 1e-8
                             mutator=1000),
               filter2dVec_in=tmpListD1,
               plotY=log10,
               plotMidX=FALSE)
d2 <-
  plotMutWaves(pdfPath=file.path(plotDir,"waves_c2_filtered_driver2_log.pdf"), 
               variants=list(drivers=0:5, superdrivers=0:7,
                             mutators=NA, passengers=NA),
               fixedParams=c(driver=0.01, factor=2, 
                             ratenormal=0.00000001, # 1e-8
                             mutator=1000),
               filter2dVec_in=tmpListD2,
               plotY=log10,
               plotMidX=FALSE)
d3 <-
  plotMutWaves(pdfPath=file.path(plotDir,"waves_c2_filtered_driver3_log.pdf"), 
               variants=list(drivers=0:5, superdrivers=0:7,
                             mutators=NA, passengers=NA),
               fixedParams=c(driver=0.01, factor=2, 
                             ratenormal=0.00000001, # 1e-8
                             mutator=1000),
               filter2dVec_in=tmpListD3,
               plotY=log10,
               plotMidX=FALSE)
d4 <-
  plotMutWaves(pdfPath=file.path(plotDir,"waves_c2_filtered_driver4_log.pdf"), 
               variants=list(drivers=0:5, superdrivers=0:7,
                             mutators=NA, passengers=NA),
               fixedParams=c(driver=0.01, factor=2, 
                             ratenormal=0.00000001, # 1e-8
                             mutator=1000),
               filter2dVec_in=tmpListD4,
               plotY=log10,
               plotMidX=FALSE)
d5 <-
  plotMutWaves(pdfPath=file.path(plotDir,"waves_c2_filtered_driver5_log.pdf"), 
               variants=list(drivers=0:5, superdrivers=0:7,
                             mutators=NA, passengers=NA),
               fixedParams=c(driver=0.01, factor=2, 
                             ratenormal=0.00000001, # 1e-8
                             mutator=1000),
               filter2dVec_in=tmpListD5,
               plotY=log10,
               plotMidX=FALSE)
#
d0_left <-
  plotMutWaves(pdfPath=file.path(plotDir,"waves_c2_filtered_driver0_log.pdf"), 
               variants=list(drivers=0:5, superdrivers=0:7,
                             mutators=NA, passengers=NA),
               fixedParams=c(driver=0.01, factor=2, 
                             ratenormal=0.00000001, # 1e-8
                             mutator=1000),
               filter2dVec_in=tmpListD0,
               plotY=log10,
               plotMidX=FALSE,
               midFunX=min)
d1_left <-
  plotMutWaves(pdfPath=file.path(plotDir,"waves_c2_filtered_driver1_log.pdf"), 
               variants=list(drivers=0:5, superdrivers=0:7,
                             mutators=NA, passengers=NA),
               fixedParams=c(driver=0.01, factor=2, 
                             ratenormal=0.00000001, # 1e-8
                             mutator=1000),
               filter2dVec_in=tmpListD1,
               plotY=log10,
               plotMidX=FALSE,
               midFunX=min)
d2_left <-
  plotMutWaves(pdfPath=file.path(plotDir,"waves_c2_filtered_driver2_log.pdf"), 
               variants=list(drivers=0:5, superdrivers=0:7,
                             mutators=NA, passengers=NA),
               fixedParams=c(driver=0.01, factor=2, 
                             ratenormal=0.00000001, # 1e-8
                             mutator=1000),
               filter2dVec_in=tmpListD2,
               plotY=log10,
               plotMidX=FALSE,
               midFunX=min)
d3_left <-
  plotMutWaves(pdfPath=file.path(plotDir,"waves_c2_filtered_driver3_log.pdf"), 
               variants=list(drivers=0:5, superdrivers=0:7,
                             mutators=NA, passengers=NA),
               fixedParams=c(driver=0.01, factor=2, 
                             ratenormal=0.00000001, # 1e-8
                             mutator=1000),
               filter2dVec_in=tmpListD3,
               plotY=log10,
               plotMidX=FALSE,
               midFunX=min)
d4_left <-
  plotMutWaves(pdfPath=file.path(plotDir,"waves_c2_filtered_driver4_log.pdf"), 
               variants=list(drivers=0:5, superdrivers=0:7,
                             mutators=NA, passengers=NA),
               fixedParams=c(driver=0.01, factor=2, 
                             ratenormal=0.00000001, # 1e-8
                             mutator=1000),
               filter2dVec_in=tmpListD4,
               plotY=log10,
               plotMidX=FALSE,
               midFunX=min)
d5_left <-
  plotMutWaves(pdfPath=file.path(plotDir,"waves_c2_filtered_driver5_log.pdf"), 
               variants=list(drivers=0:5, superdrivers=0:7,
                             mutators=NA, passengers=NA),
               fixedParams=c(driver=0.01, factor=2, 
                             ratenormal=0.00000001, # 1e-8
                             mutator=1000),
               filter2dVec_in=tmpListD5,
               plotY=log10,
               plotMidX=FALSE,
               midFunX=min)
############################
d0dist <- stepwiseDistances(d0,tmpListD0[2:6])
d1dist <- stepwiseDistances(d1,tmpListD1[1:5])
d2dist <- stepwiseDistances(d2,tmpListD2[1:5])
d3dist <- stepwiseDistances(d3,tmpListD3[1:4])
d4dist <- stepwiseDistances(d4,tmpListD4[1:4])
d5dist <- stepwiseDistances(d5,tmpListD5[1:4])

d0dist_left <- stepwiseDistances(d0_left,tmpListD0[1:8])
d1dist_left <- stepwiseDistances(d1_left,tmpListD1[1:8])
d2dist_left <- stepwiseDistances(d2_left,tmpListD2[1:8])
d3dist_left <- stepwiseDistances(d3_left,tmpListD3[1:8])
d4dist_left <- stepwiseDistances(d4_left,tmpListD4[1:8])
d5dist_left <- stepwiseDistances(d5_left,tmpListD5[1:8])

s0dist <- stepwiseDistances(s0,tmpListS0[3:7])
s1dist <- stepwiseDistances(s1,tmpListS1[1:7])
s2dist <- stepwiseDistances(s2,tmpListS2[1:7])
s3dist <- stepwiseDistances(s3,tmpListS3[1:5])
s4dist <- stepwiseDistances(s4,tmpListS3[1:3])

s0dist_left <- stepwiseDistances(s0_left,tmpListS0[2:7])
s1dist_left <- stepwiseDistances(s1_left,tmpListS1[1:7])
s2dist_left <- stepwiseDistances(s2_left,tmpListS2[1:7])
s3dist_left <- stepwiseDistances(s3_left,tmpListS3[1:8])
s4dist_left <- stepwiseDistances(s4_left,tmpListS4[1:9])
s5dist_left <- stepwiseDistances(s5_left,tmpListS5[1:9])
## rest s4, s5, etc. not enough FULL waves

#############
## boxplot ##
pdf(file.path(plotDir,"boxplot_stepwiseDistances.pdf"))
inlist <- list(d0dist,d1dist,d2dist,d3dist,d4dist,d5dist)
names(inlist) <- 0:5
boxplotMidX(inlist, main="M=1000, c=2, u=1e-8", xlab="Number of fixed drivers [r]",
            ylab="Generations between waves")
dev.off()

pdf(file.path(plotDir,"boxplot_stepwiseDistances_leftEdges.pdf"))
inlist <- list(d0dist_left,d1dist_left,d2dist_left,d3dist_left,d4dist_left,d5dist_left)
names(inlist) <- 0:5
boxplotMidX(inlist, main="M=1000, c=2, u=1e-8", xlab="Number of fixed drivers [r]",
            ylab="Generations between waves")
dev.off()

pdf(file.path(plotDir,"boxplot_stepwiseDistances_ylim.pdf"))
inlist <- list(d0dist,d1dist,d2dist,d3dist,d4dist,d5dist)
names(inlist) <- 0:5
boxplotMidX(inlist, main="M=1000, c=2, u=1e-8", xlab="Number of fixed drivers [r]",
            ylab="Generations between waves", ylim=c(0,4500))
dev.off()

## superdriver

pdf(file.path(plotDir,"boxplot_stepwiseDistances_superdriver.pdf"))
inlist <- list(s0dist,s1dist,s2dist,s3dist,s4dist)
names(inlist) <- 0:4
boxplotMidX(inlist, main="M=1000, c=2, u=1e-8", xlab="Number of fixed superdrivers [s]",
            ylab="Generations between waves")
dev.off()

pdf(file.path(plotDir,"boxplot_stepwiseDistances_superdriver_leftEdges.pdf"))
inlist <- list(s0dist_left,s1dist_left,s2dist_left,s3dist_left,s4dist_left,s5dist_left)
names(inlist) <- 0:5
boxplotMidX(inlist, main="M=1000, c=2, u=1e-8", xlab="Number of fixed superdrivers [s]",
            ylab="Generations between waves")
dev.off()

pdf(file.path(plotDir,"boxplot_stepwiseDistances_superdriver_ylim.pdf"))
inlist <- list(s0dist,s1dist,s2dist,s3dist,s4dist)
names(inlist) <- 0:4
boxplotMidX(inlist, main="M=1000, c=2, u=1e-8", xlab="Number of fixed superdrivers [s]",
            ylab="Generations between waves", ylim=c(0,4500))
dev.off()
#################################

plotMutWaves(pdfPath=file.path(plotDir,"waves_c2_filtered_driver0driver3.pdf"), 
             #              variants=list(#drivers=0:7, superdrivers=0:3,
             #                            mutators=NA, passengers=NA),
             fixedParams=c(driver=0.01, factor=2, 
                           ratenormal=0.00000001, # 1e-8
                           mutator=1000),
             filter2dVec_in=c(tmpListD0,tmpListD3))
plotMutWaves(pdfPath=file.path(plotDir,"waves_c2_filtered_driver0driver3_log.pdf"), 
             #              variants=list(#drivers=0:7, superdrivers=0:3,
             #                            mutators=NA, passengers=NA),
             fixedParams=c(driver=0.01, factor=2, 
                           ratenormal=0.00000001, # 1e-8
                           mutator=1000),
             filter2dVec_in=c(tmpListD0,tmpListD3),
             plotY=log10)


# high superdriver advantage
plotMutWaves(pdfPath=file.path(plotDir,"waves_c3.pdf"), 
             #              variants=list(#drivers=0:7, superdrivers=0:3,
             #                            mutators=NA, passengers=NA),
             fixedParams=c(driver=0.01, factor=3, 
                           ratenormal=0.00000001, # 1e-8
                           mutator=1000))
plotMutWaves(pdfPath=file.path(plotDir,"waves_c3_log.pdf"), 
             #              variants=list(#drivers=0:7, superdrivers=0:3,
             #                            mutators=NA, passengers=NA),
             fixedParams=c(driver=0.01, factor=3, 
                           ratenormal=0.00000001, # 1e-8
                           mutator=1000),
             plotY=log10)
## filtered
plotMutWaves(pdfPath=file.path(plotDir,"waves_c3_filtered_superdriver0.pdf"), 
             #              variants=list(#drivers=0:7, superdrivers=0:3,
             #                            mutators=NA, passengers=NA),
             fixedParams=c(driver=0.01, factor=3, 
                           ratenormal=0.00000001, # 1e-8
                           mutator=1000),
             filter2dVec_in=tmpListS1)
plotMutWaves(pdfPath=file.path(plotDir,"waves_c3_filtered_superdriver1.pdf"), 
             #              variants=list(#drivers=0:7, superdrivers=0:3,
             #                            mutators=NA, passengers=NA),
             fixedParams=c(driver=0.01, factor=3, 
                           ratenormal=0.00000001, # 1e-8
                           mutator=1000),
             filter2dVec_in=tmpListS4)
#
plotMutWaves(pdfPath=file.path(plotDir,"waves_c3_filtered_driver0.pdf"), 
             #              variants=list(#drivers=0:7, superdrivers=0:3,
             #                            mutators=NA, passengers=NA),
             fixedParams=c(driver=0.01, factor=3, 
                           ratenormal=0.00000001, # 1e-8
                           mutator=1000),
             filter2dVec_in=tmpListD0)
plotMutWaves(pdfPath=file.path(plotDir,"waves_c3_filtered_driver1.pdf"), 
             #              variants=list(#drivers=0:7, superdrivers=0:3,
             #                            mutators=NA, passengers=NA),
             fixedParams=c(driver=0.01, factor=3, 
                           ratenormal=0.00000001, # 1e-8
                           mutator=1000),
             filter2dVec_in=tmpListD1)
#
plotMutWaves(pdfPath=file.path(plotDir,"waves_c3_filtered_superdriver1superdriver3.pdf"), 
             #              variants=list(#drivers=0:7, superdrivers=0:3,
             #                            mutators=NA, passengers=NA),
             fixedParams=c(driver=0.01, factor=3, 
                           ratenormal=0.00000001, # 1e-8
                           mutator=1000),
             filter2dVec_in=c(tmpListS1,tmpListS3))
plotMutWaves(pdfPath=file.path(plotDir,"waves_c3_filtered_superdriver1superdriver3_log.pdf"), 
             #              variants=list(#drivers=0:7, superdrivers=0:3,
             #                            mutators=NA, passengers=NA),
             fixedParams=c(driver=0.01, factor=3, 
                           ratenormal=0.00000001, # 1e-8
                           mutator=1000),
             filter2dVec_in=c(tmpListS1,tmpListS3),
             plotY=log10)
plotMutWaves(pdfPath=file.path(plotDir,"waves_c3_filtered_driver0driver3.pdf"), 
             #              variants=list(#drivers=0:7, superdrivers=0:3,
             #                            mutators=NA, passengers=NA),
             fixedParams=c(driver=0.01, factor=3, 
                           ratenormal=0.00000001, # 1e-8
                           mutator=1000),
             filter2dVec_in=c(tmpListD0,tmpListD3))
plotMutWaves(pdfPath=file.path(plotDir,"waves_c3_filtered_driver0driver3_log.pdf"), 
             #              variants=list(#drivers=0:7, superdrivers=0:3,
             #                            mutators=NA, passengers=NA),
             fixedParams=c(driver=0.01, factor=3, 
                           ratenormal=0.00000001, # 1e-8
                           mutator=1000),
             filter2dVec_in=c(tmpListD0,tmpListD3),
             plotY=log10)
# 
# # no superdriver advantage
# plotMutWaves(pdfPath=file.path(plotDir,"waves_c1.pdf"), 
#              #              variants=list(#drivers=0:7, superdrivers=0:3,
#              #                            mutators=NA, passengers=NA),
#              fixedParams=c(driver=0.01, factor=1, 
#                            ratenormal=0.00000001, # 1e-8
#                            mutator=1000))
# plotMutWaves(pdfPath=file.path(plotDir,"waves_c1_log.pdf"), 
#              #              variants=list(#drivers=0:7, superdrivers=0:3,
#              #                            mutators=NA, passengers=NA),
#              fixedParams=c(driver=0.01, factor=1, 
#                            ratenormal=0.00000001, # 1e-8
#                            mutator=1000),
#              plotY=log10)
# plotMutWaves(pdfPath=file.path(plotDir,"waves_c1_filtered_superdriver0.pdf"), 
#              #              variants=list(#drivers=0:7, superdrivers=0:3,
#              #                            mutators=NA, passengers=NA),
#              fixedParams=c(driver=0.01, factor=1, 
#                            ratenormal=0.00000001, # 1e-8
#                            mutator=1000),
#              filter2dVec_in=tmpListS1)
# plotMutWaves(pdfPath=file.path(plotDir,"waves_c1_filtered_superdriver1.pdf"), 
#              #              variants=list(#drivers=0:7, superdrivers=0:3,
#              #                            mutators=NA, passengers=NA),
#              fixedParams=c(driver=0.01, factor=1, 
#                            ratenormal=0.00000001, # 1e-8
#                            mutator=1000),
#              filter2dVec_in=tmpListS4)
# #
# plotMutWaves(pdfPath=file.path(plotDir,"waves_c1_filtered_driver0.pdf"), 
#              #              variants=list(#drivers=0:7, superdrivers=0:3,
#              #                            mutators=NA, passengers=NA),
#              fixedParams=c(driver=0.01, factor=1, 
#                            ratenormal=0.00000001, # 1e-8
#                            mutator=1000),
#              filter2dVec_in=tmpListD0)
# plotMutWaves(pdfPath=file.path(plotDir,"waves_c1_filtered_driver1.pdf"), 
#              #              variants=list(#drivers=0:7, superdrivers=0:3,
#              #                            mutators=NA, passengers=NA),
#              fixedParams=c(driver=0.01, factor=1, 
#                            ratenormal=0.00000001, # 1e-8
#                            mutator=1000),
#              filter2dVec_in=tmpListD1)
