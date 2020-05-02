

readStoreCasePop <- function(case, parallel=goParallel, readPopParallel=goParallel) {
  namesOfCase <- grep(case, popsFileNames, value=T)
  fullNamesOfCase <- sapply(namesOfCase, function(x) file.path(populationsDir, x))
#   for (dir in fullNamesOfCase) {
#     basenameDir <- basename(dir)
#     tmpName <- file.path(rdataDir,sprintf("%s.rda",basenameDir))
#     if (!file.exists(tmpName)) {
#       maxRepsExists <- file.exists(file.path(dir,sprintf("r%s.pop",maxReps))) # e.g. r050.pop
#       if (maxRepsExists) { # if dir not yet finished (i.e. maxReps not in dir)
#         print(dir)
#         tmpPop <- readPopOfDir(dir, "t.+.pop", parallel=parallel, cores=nCores)
#         if (!is.null(tmpPop)) { # directory was not empty (comment: really needed if checked for r050?)
#           #       tmpPop <- readPopOfDir(dir, "t00[1-3]+-g04491.pop", parallel=T)
#           tmpPop <- separateIntermediatePops(tmpPop, 'g')
#           basenameDir_trunc <- basenameDir #sub(".+\\.rda", "\\1", basenameDir)
#           assign(basenameDir_trunc, tmpPop)
#           #       tmpObject <- get(basenameDir)
#           save(list=c(basenameDir_trunc), file=tmpName)
#           rm(basenameDir, tmpPop, basenameDir_trunc)
#         }
#       }       
#     } else {
#       warning(sprintf("populations %s already saved to rdata. ignoring directory..",dir))
#       #       load(tmpName)
#     }
#   }
  require(parallel)
   mclapply(fullNamesOfCase, function(dir) {
    basenameDir <- basename(dir)
    tmpName <- file.path(rdataDir,sprintf("%s.rda",basenameDir))
    if (!file.exists(tmpName)) {
      maxRepsExists <- file.exists(file.path(dir,sprintf("r%s.pop",maxReps))) # e.g. r050.pop
      if (maxRepsExists) { # if dir not yet finished (i.e. maxReps not in dir)
        print(dir)
#         readPopParallel <- ifelse(grepl("mutator1",dir), FALSE, TRUE)
        if (grepl("mutator1",dir))
          readPopParallel <- FALSE
        tmpPop <- readPopOfDir(dir, "t.+.pop", parallel=readPopParallel, cores=nCores)
        if (!is.null(tmpPop)) { # directory was not empty (comment: really needed if checked for r050?)
          #       tmpPop <- readPopOfDir(dir, "t00[1-3]+-g04491.pop", parallel=T)
          tmpPop <- separateIntermediatePops(tmpPop, 'g')
          basenameDir_trunc <- basenameDir #sub(".+\\.rda", "\\1", basenameDir)
          assign(basenameDir_trunc, tmpPop)
          #       tmpObject <- get(basenameDir)
          save(list=c(basenameDir_trunc), file=tmpName)
          rm(basenameDir, tmpPop, basenameDir_trunc)
        }
      }       
    } else {
      warning(sprintf("populations %s already saved to rdata. ignoring directory..",dir))
      #       load(tmpName)
    }
  }, mc.cores=nCores)
}

readStoreCasePop_finalPops <- function(case, parallel=goParallel) {
  namesOfCase <- grep(case, popsFileNames, value=T)
  fullNamesOfCase <- sapply(namesOfCase, function(x) file.path(populationsDir, x))
  print(fullNamesOfCase)
  for (dir in fullNamesOfCase) {
    basenameDir <- basename(dir)
    tmpName <- file.path(rdataDir,sprintf("finalPops_%s.rda",basenameDir))
    print(tmpName)
    #     basenameDir <- basename(basenameDir)
    if (!file.exists(tmpName)) {
      maxRepsExists <- file.exists(file.path(dir,sprintf("r%s.pop",maxReps)))
      if (maxRepsExists) { # if dir not yet finished (i.e. maxReps not in dir)
        print(dir)
        tmpPop <- readPopOfDir(dir, "r.+.pop", parallel=parallel, cores=nCores)
        if (!is.null(tmpPop)) { # directory was not empty 
          #       tmpPop <- readPopOfDir(dir, "t00[1-3]+-g04491.pop", parallel=T)
          basenameDir_trunc <- sprintf("finalPops_%s",basenameDir) #sub(".+\\.rda", "\\1", basenameDir)
          assign(basenameDir_trunc, tmpPop)
          #       tmpObject <- get(basenameDir)
          save(list=c(basenameDir_trunc), file=tmpName)
          rm(basenameDir, tmpPop, basenameDir_trunc)
        }
      }       
    } else {
      warning(sprintf("populations %s already saved to rdata. ignoring directory..",dir))
      #       load(tmpName)
    }
  }
}



#   pops_local <- 
#     readPopOfDirs(, "t.+.pop", parallel=F)
# #   pops_separated_local <- separateIntermediatePops2(pops_local, 'g')
#   pops_local

# tmp_A_name <- file.path(rdataDir,"popsA.rda")
# if (!file.exists(tmp_A_name)) {
#   popsA <- readCasePop("A")
#   print("read in case A")
#   save(popsA, file=tmp_A_name)
#   print("intermediate save")
# } else {
#   warning("populations A already save in rdata. loading file..")
#   load(tmp_A_name)
# }
# 
# # stop("debug")
#    
# tmp_B_name <- file.path(rdataDir,"popsB.rda")
# if (!file.exists(tmp_B_name)) {
#   popsB <- readCasePop("B")
#   print("read in case B")
#   save(popsB, file=tmp_B_name)
#   print("intermediate save")
# } else  {
#   warning("populations B already save in rdata. loading file..")
#   load(tmp_B_name)
# }
# 
# tmp_C_name <- file.path(rdataDir,"popsC.rda")
# if (!file.exists(tmp_C_name)) {
#   popsC <- readCasePop("C")
#   print("read in case C")
#   save(popsC, file=tmp_C_name)
#   print("intermediate save")
# } else  {
#   warning("populations C already save in rdata. loading file..")
#   load(tmp_C_name)
# }
# 
# tmp_D_name <- file.path(rdataDir,"popsD.rda")
# if (!file.exists(tmp_D_name)) {
#   popsD <- readCasePop("D")
#   print("read in case D")
#   save(popsD, file=tmp_D_name)
#   print("intermediate save")
# } else  {
#   warning("populations D already save in rdata. loading file..")
#   load(tmp_D_name)
# }
# 
# #   pops <- readPopOfDirs(sapply(popsFileNames, function(x) file.path(populationsDir, x)), "t.+.pop")
# print("separating pops..")
# popsA_separated <- separateIntermediatePops2(popsA, 'g')
# print("separated A")
# popsB_separated <- separateIntermediatePops2(popsB, 'g')
# print("separated B")
# popsC_separated <- separateIntermediatePops2(popsC, 'g')
# print("separated C")
# popsD_separated <- separateIntermediatePops2(popsD, 'g')
# print("separated D")
# 
# print("saving populations to RData")
# save(c("popsA_separated",
#        "popsB_separated",
#        "popsC_separated",
#        "popsD_separated"), 
#      file=pops_RDataName)
# print("removing initial pops object")
# rm(popsA, popsB, popsC, popsD)


