############################################################################
### this file contains functions to pre-analyse finished simulation results
### (i.e., contains a corresponding rdata file).
###
### Efficiency of analysis is increased by storing intermediate results and 
### loading them if existing. Functions can be re-used to update with newly
### finished simulation results (i.e., new rdata file in rdataDir)-
###
### first it is checked, whether or not a corresponding rdata file exists
### in rdataDir, then it is loaded (if not yet in workspace!), then it is checked
### whether or not full analysis has been performed and stored to intermediateDir already.
### If so, rdata is loaded. if not, analysis steps are performed and stored. intermediate 
### results are checked in the same way (if existing, not performed, but loaded)
### steps are performed and intermediate results are stored.
### (see getAnalysis.generic for initiation and further details)
############################################################################

#' use function to perform an analysis step only if result not stored
#' @param obj usually a list of list of tabs corresponding to list of
#' generations of replicates of populations.
#' @param objName name that result analysis on obj should get in environment
#' @param fun function to apply to obj
#' @param ... additional parameters to fun
#' @return result of fun(obj) named as objName
performAnalysis <- function(obj, objName, fun, ...) {
  fn_path <- file.path(intermediateDir, sprintf("%s.rda",objName))
  if (!file.exists(fn_path)) {
    tmpRes <- fun(obj, ...)
    assign(objName, tmpRes)
    save(list=c(objName), file=fn_path)
    return(tmpRes)
  } else {
    load(fn_path)
    return(get(objName))
  }
}

#' get different types of "one clones" for each replicate in simulation result
#' full result stored in realName_oneClonesRes.rda
#' @param realName name of simulation result
#' @param ... additional arguments to innerF (most likely not used)
#' @return list of elements wav, wsd, highest, corresponding to different types
#' of "one clones". each element stores generations and replicates
forOne_oneClones <- function(realName, ...) {
  innerF <- function(tmpObject) {
    tmpOneClones_wav_name <- sprintf("%s_oneClones_wav", realName)
    tmpOneClones_wsd_name <- sprintf("%s_oneClones_wsd", realName)
    tmpOneClones_highest_name <- sprintf("%s_oneClones_highest", realName)
    
    print(sprintf("perform oneClones wav on %s", realName))
    resWav <- performAnalysis(tmpObject, tmpOneClones_wav_name, get1Clone2, which="wav") 
    
    print(sprintf("perform oneClones wsd on %s", realName))
    resWsd <- performAnalysis(tmpObject, tmpOneClones_wsd_name, get1Clone2, which="wsd") 
    
    print(sprintf("perform oneClones highest on %s", realName))
    resHighest <- performAnalysis(tmpObject, tmpOneClones_highest_name, get1Clone2, which="highestFreq") 
    
    list(wav=resWav,
         wsd=resWsd,
         highest=resHighest)
  }
  forOne(realName, "oneClonesRes", innerF, ...)
}

#' perform oneClone analysis on only final population
forOne_oneClones_finalPops <- function(realName, ...) {
  innerF <- function(tmpObject) {
    tmpOneClones_wav_name <- sprintf("%s_oneClonesFinalPops_wav", realName)
    tmpOneClones_wsd_name <- sprintf("%s_oneClonesFinalPops_wsd", realName)
    tmpOneClones_highest_name <- sprintf("%s_oneClonesFinalPops_highest", realName)
    
    print(sprintf("perform oneClones_finalPops wav on %s", realName))
    resWav <- performAnalysis(tmpObject, tmpOneClones_wav_name, get1Clone1, which="wav") 
    
    print(sprintf("perform oneClonesFinalPops wsd on %s", realName))
    resWsd <- performAnalysis(tmpObject, tmpOneClones_wsd_name, get1Clone1, which="wsd") 
    
    print(sprintf("perform oneClonesFinalPops highest on %s", realName))
    resHighest <- performAnalysis(tmpObject, tmpOneClones_highest_name, get1Clone1, which="highestFreq") 
    
    list(wav=resWav,
         wsd=resWsd,
         highest=resHighest)
  }
  forOne(realName, "oneClonesFinalPopsRes", innerF, ...)
}

#' get 4mutClasses for all generations and replicates. see comments above
forOne_4mutClasses <- function(realName, ...) {
  innerF <- function(tmpObject) {
    tmp4MutClasses_name <- sprintf("%s_4mutClasses", realName)
    N4MutClasses_name <- sprintf("%s_N4mutClasses", realName)
    print(sprintf("performing 4mutClasses on %s", realName))
    fourMutClasses <- performAnalysis(tmpObject, tmp4MutClasses_name, get4MutationClasses2,
                                      weighted=TRUE, cores=1, ...)     
    print(sprintf("performing N4mutClasses on %s", realName))
    N4mutClasses <- performAnalysis(fourMutClasses, N4MutClasses_name, getNumber4mutClasses2)    
    list(fourMutClasses=fourMutClasses,
         N4mutClasses=N4mutClasses)
  }
  forOne(realName, "4mutClassesRes", innerF, ...)
}

#' get means of 4mutClasses for each generation. see comments above
forOne_4mutClassesMeans <- function(realName, ...) {
  innerF <- function(tmpObject) {
    tmp4MutClasses_name <- sprintf("%s_4mutClassesMean", realName)
    N4MutClasses_name <- sprintf("%s_N4mutClassesMean", realName)
    print(sprintf("performing 4mutClasses mean on %s", realName))    
    fourMutClasses <- performAnalysis(tmpObject, tmp4MutClasses_name, get4MutationClasses_mean2,
                    weighted=TRUE, cores=3, ...)     
    print(sprintf("performing N4mutClasses mean on %s", realName))
    N4mutClasses <- performAnalysis(fourMutClasses, N4MutClasses_name, getNumber4mutClasses2)  
    list(fourMutClasses=fourMutClasses,
         N4mutClasses=N4mutClasses)
  }
  forOne(realName, "4mutClassesMeansRes", innerF, ...)
}

#' for final pops
forOne_4mutClassesMeans_finalPops <- function(realName, ...) {
  innerF <- function(tmpObject) {
    tmp4MutClasses_name <- sprintf("%s_4mutClassesMean_finalPops", realName)
    N4MutClasses_name <- sprintf("%s_N4mutClassesMean_finalPops", realName)
    print(sprintf("performing 4mutClasses mean on final Pops %s", realName))
    fourMutClasses <- performAnalysis(tmpObject, tmp4MutClasses_name, get4MutationClasses_mean1,
                    weighted=TRUE)    
    print(sprintf("performing N4mutClasses mean on final Pops %s", realName))
    N4mutClasses <- performAnalysis(fourMutClasses, N4MutClasses_name, getNumber4mutClasses1)  
    list(fourMutClasses=fourMutClasses,
         N4mutClasses=N4mutClasses)
  }
  forOne(realName, "4mutClassesMeanResFinalPops", innerF, ...)
}

#' for final pops
forOne_4mutClasses_finalPops <- function(realName, ...) {
  innerF <- function(tmpObject) {
    tmp4MutClasses_name <- sprintf("%s_4mutClasses_finalPops", realName)
    N4MutClasses_name <- sprintf("%s_N4mutClasses_finalPops", realName)
    print(sprintf("performing 4mutClasses on final Pops %s", realName))
    fourMutClasses <- performAnalysis(tmpObject, tmp4MutClasses_name, get4MutationClasses1,
                    weighted=TRUE, ...)  
    print(sprintf("performing N4mutClasses on final Pops %s", realName))
    N4mutClasses <- performAnalysis(fourMutClasses, N4MutClasses_name, getNumber4mutClasses1)  
    list(fourMutClasses=fourMutClasses,
         N4mutClasses=N4mutClasses)
  }
  forOne(realName, "4mutClassesResFinalPops", innerF, ...)
}

#' get waiting times for simulation results.
#' will compute waiting time and probabilites to cancer and to mutator
#' phenotype. each of this four intermediate results are stored separately.
#' also computes waiting times on switched populations, i.e. list of
#' replicates of generations of populations.
#' @param realName character name of simu result
#' @param kCancer integer for when cancer is defined
#' @return list of waiting time results
forOne_waitingTime <- function(realName, kCancer=k_cancer) {
  
  innerF <- function(tmpObject) {
    ### analysis steps ###
    
    tmpCancerWaiting <- sprintf("%s_t_cancer%s",realName,kCancer)
    tmpMutatorWaiting <- sprintf("%s_t_mutator",realName)
    tmpCancerWaiting2 <- sprintf("%s_t_cancer%s_switched",realName,kCancer)
    tmpMutatorWaiting2 <- sprintf("%s_t_mutator_switched",realName)
    
    tmpCancerWaitingProbs <- sprintf("%s_t_cancer%s_probs",realName,kCancer)
    tmpMutatorWaitingProbs <- sprintf("%s_t_mutator_probs",realName)
    
    ### waiting time to cancer ###
    
    print(sprintf("performing cancer waiting time to %s", realName))
    resCancerWaiting <- performAnalysis(tmpObject, tmpCancerWaiting, getWaitingTime_toCancer, k=kCancer) 
    
    print(sprintf("performing cancer waiting time switched to %s", realName))
    resCancerWaiting2 <- performAnalysis(switchPops2(tmpObject), tmpCancerWaiting2, getWaitingTime_toCancer, k=kCancer) 
    
    print(sprintf("performing cancer waiting time probs to %s", realName))
    resCancerWaitingProbs <- 
      performAnalysis(tmpObject, tmpCancerWaitingProbs, getWaitingTimesProb_toCancer2, k=kCancer)  
    
    ### waiting time to mutator ###
    
    print(sprintf("performing mutator waiting time to %s", realName))
    M <- as.numeric(sub(".+mutator([0-9]+).+", "\\1", realName))
    resMutatorWaiting <- performAnalysis(tmpObject, tmpMutatorWaiting, getWaitingTime_toMutator, M=M)  
    resMutatorWaiting2 <- performAnalysis(switchPops2(tmpObject), tmpMutatorWaiting2, getWaitingTime_toMutator, M=M)  
    
    print(sprintf("performing mutator waiting time probs to %s", realName))
    resMutatorWaitingProbs <- 
      performAnalysis(tmpObject, tmpMutatorWaitingProbs, getWaitingTimesProb_toMutator2, M=M)
    
    ### return results ###
    
    res <- list(cancerWaitingTime=resCancerWaiting,
                cancerWaitingTimeSwitched=resCancerWaiting2,
                mutatorWaitingTime=resMutatorWaiting,
                mutatorWaitingTimeSwitched=resMutatorWaiting2,
                cancerWaitingTimeProbs=resCancerWaitingProbs,
                mutatorWaitingTimeProbs=resMutatorWaitingProbs)
    return(res)
  }
  
  forOne(realName, "waitingTimeRes", innerF)
}

#' generic function to apply on a single simulation result
#' rdata file is loaded from rdataDir according to realName
#' this function applies innerFun on a single simu result.
#' if finalPops object is smaller than 50, it is truncated
#' to 20, if normal (inter) object's first element is smaller than
#' 50 it is truncated to 20.
#' @param realName character name of simulation result
#' @param resName character name that should be assigned to object
#' result will also be named resName.rda
#' @param innerFun function to be applied to one simulation result
#' @return result of innerFun on loaded realName
forOne <- function(realName, resName, innerFun) {
  
  tmpRealName_finalRes_name <- sprintf("%s_%s", realName, resName)
  tmpRealName_finalRes_fn <- file.path(intermediateDir, sprintf("%s.rda", tmpRealName_finalRes_name))
  
  if (!file.exists(tmpRealName_finalRes_fn)) {
    print(sprintf("current: %s is computed", tmpRealName_finalRes_fn))
    #       load(file.path(rdataDir,sprintf("%s.rda",realName)))
    if (loadAllDataGlobal) env <- .GlobalEnv  else env <- environment()
    #       env <- ifelse(loadAllDataGlobal, .GlobalEnv, environment())
    # load corresponding rdata file into workspace if not loaded already!
    print(file.path(rdataDir,sprintf("%s.rda",realName)))
    tmpObject <- loadIfNot(file.path(rdataDir,sprintf("%s.rda",realName)), env=env)
#     tmpObject <- get(realName, envir=env)
#     if (grepl("^finalPops",realName)) {
#       if (length(tmpObject) < desiredReps)
#           tmpObject <- truncateObjectFinal(tmpObject, minReps)      
#     }
# 
#     if (grepl("^[ABCD]_",realName)) {
#       if (length(tmpObject) > 0 && length(tmpObject))
#     }
    res <- innerFun(tmpObject)
    assign(tmpRealName_finalRes_name, res)
    print(tmpRealName_finalRes_name)
    print(tmpRealName_finalRes_fn)
    save(list=c(tmpRealName_finalRes_name), file=tmpRealName_finalRes_fn)
    return(res)
  } else {
    load(tmpRealName_finalRes_fn)
    return(get(tmpRealName_finalRes_name))
  }
}

#' generic analysis function. performs analysis on a set of simulation results corresponding
#' to a case. analysis can be performed for each simulation result in parallel (max <cores> threads)
#' @param character case specified as regex pattern which specifies simulation name
#' @param thisFun analysis function to be applied on each simulation result
#' @param ... additional arguments to thisFun
#' @param parallel logical if analysis should be parallized
#' @param cores integer specifies max number of cores to be used in parallel version
getAnalysis.generic <- function(case, thisFun, ..., parallel=goParallel, cores=nCores) {
  rdaNames <- grep(sprintf("%s.*\\.rda",case), list.files(file.path(rdataDir)), value=T)
  realNames <- sub("(.+)\\.rda", "\\1", rdaNames)
  if (!parallel) { 
    tmpList <- lapply(realNames, thisFun, ...)
  } #sapply(rdaNames, load)
  else {
    require(parallel)
    tmpList <- mclapply(realNames, thisFun, ..., mc.cores=nCores, mc.preschedule=FALSE)
  } 
  names(tmpList) <- realNames
  tmpList
}

getAnalysis.waitingTimes <- function(case, parallel=goParallel, cores=nCores, kCancer=k_cancer)
  getAnalysis.generic(case=case, parallel=parallel, cores=cores, thisFun=forOne_waitingTime, kCancer=kCancer)

getAnalysis.4mutClassesMeans <- function(case, parallel=goParallel, cores=nCores)
  getAnalysis.generic(case=case, parallel=parallel, cores=cores, thisFun=forOne_4mutClassesMeans)

getAnalysis.4mutClassesMeans_finalPops <- function(case, parallel=goParallel, cores=nCores)
  getAnalysis.generic(case=case, parallel=parallel, cores=cores, thisFun=forOne_4mutClassesMeans_finalPops)

getAnalysis.4mutClasses <- function(case, parallel=goParallel, cores=nCores)
  getAnalysis.generic(case=case, parallel=parallel, cores=cores, thisFun=forOne_4mutClasses)

getAnalysis.4mutClasses_finalPops <- function(case, parallel=goParallel, cores=nCores)
  getAnalysis.generic(case=case, parallel=parallel, cores=cores, thisFun=forOne_4mutClasses_finalPops)

getAnalysis.oneClones <-function(case, parallel=goParallel, cores=nCores)
  getAnalysis.generic(case=case, parallel=parallel, cores=cores, thisFun=forOne_oneClones)

getAnalysis.oneClones_finalPops <-function(case, parallel=goParallel, cores=nCores)
  getAnalysis.generic(case=case, parallel=parallel, cores=cores, thisFun=forOne_oneClones_finalPops)
