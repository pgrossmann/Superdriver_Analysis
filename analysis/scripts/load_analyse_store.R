loadIfNot <- function(name, env=.GlobalEnv) {
  if (! name%in%ls(env))
    load(file.path(rdataDir,sprintf("%s.rda",name)))
  assign(name, get(name), env)
  return(get(name))
}

getAnalysis_waitingTimes <- function(case, parallel=goParallel, cores=nCores, kCancer=k_cancer) {
  rdaNames <- grep(sprintf("%s.*\\.rda",case), list.files(file.path(rdataDir)), value=T)
  realNames <- sub("(.+)\\.rda", "\\1", rdaNames)
  intermediateDir <- file.path(rdataDir, "intermediate")
  if (!file.exists(intermediateDir))
    dir.create(intermediateDir, showWarnings = FALSE)
  
  forOne <- function(realName) {
    print(sprintf("current: %s", realName))
    
    tmpRealName_finalRes_name <- sprintf("%s_waitingTimeRes",realName)
    tmpRealName_finalRes_fn <- file.path(intermediateDir, sprintf("%s.rda", tmpRealName_finalRes_name))
    
    if (!file.exists(tmpRealName_finalRes_fn)) {
#       load(file.path(rdataDir,sprintf("%s.rda",realName)))
      if (loadAllDataGlobal) env <- .GlobalEnv  else env <- environment()
#       env <- ifelse(loadAllDataGlobal, .GlobalEnv, environment())
      tmpObject <- loadIfNot(realName, env=env) # load corresponding rdata file into workspace if not loaded already!
      tmpObject <- get(realName, envir=env)
      ### analysis steps ###
      
      tmpCancerWaiting <- sprintf("%s_t_cancer",realName)
      tmpMutatorWaiting <- sprintf("%s_t_mutator",realName)
      
      tmpCancerWaitingProbs <- sprintf("%s_t_cancer_probs",realName)
      tmpMutatorWaitingProbs <- sprintf("%s_t_mutator_probs",realName)
      
      #' use function to perform an analysis step only if result not stored
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
      
      ### waiting time to cancer ###
      
      print(sprintf("performing cancer waiting time to %s", realName))
      resCancerWaiting <- performAnalysis(tmpObject, tmpCancerWaiting, getWaitingTime_toCancer, k=kCancer)  
      
      print(sprintf("performing cancer waiting time probs to %s", realName))
      resCancerWaitingProbs <- 
        performAnalysis(tmpObject, tmpCancerWaitingProbs, getWaitingTimesProb_toCancer2, k=kCancer)  
      
      ### waiting time to mutator ###
      
      print(sprintf("performing mutator waiting time to %s", realName))
      M <- as.numeric(sub(".+mutator([0-9]+).+", "\\1", realName))
      resMutatorWaiting <- performAnalysis(tmpObject, tmpMutatorWaiting, getWaitingTime_toMutator, M=M)  
      
      print(sprintf("performing mutator waiting time probs to %s", realName))
      resMutatorWaitingProbs <- 
        performAnalysis(tmpObject, tmpMutatorWaitingProbs, getWaitingTimesProb_toMutator2, M=M)
      
      ### return results ###
      
      res <- list(cancerWaitingTime=resCancerWaiting,
                  mutatorWaitingTime=resMutatorWaiting,
                  cancerWaitingTimeProbs=resCancerWaitingProbs,
                  mutatorWaitingTimeProbs=resMutatorWaitingProbs)
      assign(tmpRealName_finalRes_name, res)
      save(list=c(tmpRealName_finalRes_name), file=tmpRealName_finalRes_fn)
      return(res)
    } else {
      load(tmpRealName_finalRes_fn)
      return(get(tmpRealName_finalRes_name))
    }
  }
  
  if (!parallel) {
    tmpList <- lapply(realNames, forOne)
  } #sapply(rdaNames, load)
  else {
    require(parallel)
    tmpList <- mclapply(realNames, forOne, mc.cores=nCores, mc.preschedule=FALSE)
  } 
  names(tmpList) <- realNames
  tmpList
}


