print("loading rdata..")

#' case character specifying which case (A-D) to load
loadAndTransform <- function(case, parallel=goParallel) {
  
  ## read in ##
  rdaNames <- grep(sprintf("%s.+\\.rda",case), list.files(file.path(rdataDir)), value=T)
  realNames <- sub("(.+)\\.rda", "\\1", rdaNames)
  if (!parallel) {
    tmpList <- lapply(realNames, function(realName) {
      load(file.path(rdataDir,sprintf("%s.rda",realName)))
      get(realName)
    })
  } #sapply(rdaNames, load)
  else {
    require(parallel)
    tmpList <- mclapply(realNames, function(realName) {
      load(file.path(rdataDir,sprintf("%s.rda",realName)))
      get(realName)
    }, mc.cores=nCores, mc.preschedule=FALSE)
  } 
#   tmpList <- lapply(realNames, get) # object names in each rda file equals name of file without.rda
  names(tmpList) <- realNames
  
  ## chose parts of data ##
  #########################
  
  ## return ##
  
#   assign(sprintf("%s_pops", case), tmpList)
  tmpList
}


