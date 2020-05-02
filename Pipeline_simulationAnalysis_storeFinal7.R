##################################################
###
### 
### patrick@jimmy.harvard.edu
###
### execute within ./R/ directory!
##################################################

#--------------------
### prerequisites ###
#--------------------

rm(list=ls(all=TRUE))

#-------------------
### requirements ###
#-------------------



#-----------------------
### global variables ###
#-----------------------

scriptDir <- file.path("analysis","scripts")
rdataDir <- file.path("/links/grid/scratch/grossman","rdata")
plotDir <- file.path("analysis","plots")
sourcef <- function(file) source(file.path(scriptDir, file))
toPlots <- function(file) file.path(plotDir,file)

populationsDir <- 
  "/links/grid/scratch/grossman/simulation_out/finalRuns/4500gens_N10e9_n10e6ratenormal1e-7"
popsFileNames <- list.files(populationsDir)
# file.path("~","simulation_out","finalRuns","4500gens_N10e9_n10e6normalrate1e-8") 
#"~/simulation_out/finalRuns/4500gens_N10e9_n10e6" # file.path()
# file.path("~","simulation_out","finalRuns","4500gens_N10e9_n10e6")

intermediateDir <- file.path(rdataDir, "intermediate")
if (!file.exists(intermediateDir))
  dir.create(intermediateDir, showWarnings = FALSE)

finalPopsDir <- file.path(rdataDir, "finalPops")
if (!file.exists(finalPopsDir))
  dir.create(finalPopsDir, showWarnings = FALSE)

finalPopsDir <- intermediateDir

rdataFiles <- c()
nCores <- 50
goParallel <- TRUE
infMutators <- 1000 # constant used in clonex!
maxReps <- "020" # constant used in clonex!

testYes <- FALSE
storePopsYes <- FALSE
storeFinalPopsYes <- FALSE
readPopsYes <- FALSE
preAnalyseYes <- FALSE
analyseMoreYes <- FALSE
plotYes <- FALSE
loadAllData <- FALSE # used for calls to loadIfNot()

#  testYes <- TRUE
#  storePopsYes <- TRUE
  storeFinalPopsYes <- TRUE
#  readPopsYes <- TRUE
#  preAnalyseYes <- TRUE
#  analyseMoreYes <- TRUE
#  plotYes <- TRUE
loadAllDataGlobal <- TRUE 

k_cancer <- 20

#----------------
### load data ###
#----------------



#----------------------------
### load common functions ###
#----------------------------

wd_tmp <- getwd()
# setwd("~/Dropbox/current_projects/717genesignatures/R/analysis/scripts")
# source("custom_functions_only.R")
setwd(wd_tmp)

##

sourcef("analysis_functions.R")
sourcef("plotting_functions.R")

#--------------------------------------
### additional function definitions ###
#--------------------------------------

#' @param name string full path to rda file including .rda
loadIfNot <- function(name, env=.GlobalEnv) {
  nameLoad <- name
  name <- sub("(.+)\\.rda", "\\1", basename(name))
  if (! name%in%ls(env))
    load(nameLoad)
  else
    print("not loading")
  assign(name, get(name), env)
  return(get(name))
}

#' get percentage of u most occurend genes within a set of gene signatures, input is a dendrogram x
#' @param g list of gene signatures
#' @param u first u most occurend genes
#' @param d dendrogram to obtain genes from (dendro is meant to be sub dendro of GeneSigDB signatures)
#' @return allGenes list of genes and their percentage of occurence within the dendrogram cluster

#######################################################################################################################################

print("prerequisite done..")

#################################################################################
######################### ##### prepare analysis  ##### #########################
#################################################################################

if (testYes) {
  
  sourcef("testAnalysisCode.R")
  stop("only testing so far")
  rm(list=ls(all=TRUE))
}


#------------------------
### store populations ###
#------------------------
if (storePopsYes) {
  
  pops_RDataName <- file.path(rdataDir,"pops_separated.rda")
  rdataFiles <- c(rdataFiles, pops_RDataName)
  
  if (!file.exists(pops_RDataName)) {
    sourcef("storePopulations.R")
    system.time(sourcef("storePopulationsExecute.R"))
  } else load(pops_RDataName)
  
  rm(list=ls(all=TRUE))
  stop("stored intermediate pops")
  
}

if (storeFinalPopsYes) {
  sourcef("storePopulations.R")
  #   readStoreCasePop_finalPops("A_mutator1000driver0.001factor1")
  #   stop("testest")
  print("reading final pops A")
  readStoreCasePop_finalPops("A")
  print("reading final pops B")
  readStoreCasePop_finalPops("B")
  print("reading final pops C")
  readStoreCasePop_finalPops("C")
  print("reading final pops D")
  readStoreCasePop_finalPops("D")
  stop("stored final pops")
}

#-----------------------
### load populations ###
#-----------------------

## not used anymore!!
if (readPopsYes) {
  
  sourcef("readPopulations.R")
  
  popsA <- loadAndTransform("A", parallel=FALSE)
  print("case A is ready!")
  popsB <- loadAndTransform("B", parallel=FALSE)
  print("case B is ready!")
  popsC <- loadAndTransform("C", parallel=FALSE)
  print("case C is ready!")
  popsD <- loadAndTransform("D", parallel=FALSE)
  print("case D is ready!")
  
  print("Population data ready!")
  
}


###############################################################################
######################### ##### start analysis  ##### #########################
###############################################################################

####################### ##### pre-analyse populations  ##### ##################
### this block performs initial analyses steps, such as waiting time, mutation classes,
### or computation of weighted averages clones on already finished populations and 
### stores intermediate results in order to avoid re-computation

if (preAnalyseYes) {
  sourcef("pre-analyseANDstore.R")
} 

#################### ##### interpret pre-analyses results  ##### ###############
### this block analyses and interprets finished pre-analysis results, plots figures
### to files, and perhaps stores further results

if (analyseMoreYes) {
  sourcef("analyseANDplot.R")
} else stop("Pipeline end")

##########################
### waiting time analysis
##########################



##############
### Diversity
##############

#---------------------
### final populations 
#---------------------

#----------------------------
### intermediate populations 
#----------------------------





####################################################################

stop("below testing area")

##########################
### waiting time analysis
##########################

#-----------------------------
### waiting time to cancer ###
#-----------------------------

waitingTimes_cancer_popsA <- getWaitingTimes_toCancer_store(popsA, k_cancer, store=TRUE)
waitingTimes_cancer_popsB <- getWaitingTimes_toCancer_store(popsB, k_cancer, store=TRUE)
waitingTimes_cancer_popsC <- getWaitingTimes_toCancer_store(popsC, k_cancer, store=TRUE)
waitingTimes_cancer_popsD <- getWaitingTimes_toCancer_store(popsD, k_cancer, store=TRUE)

### probabilities

waitingTimes_cancer_probs_popsA <- getWaitingTimesProb_toCancer3_store(popsA, k_cancer, store=TRUE)
waitingTimes_cancer_probs_popsB <- getWaitingTimesProb_toCancer3_store(popsB, k_cancer, store=TRUE)
waitingTimes_cancer_probs_popsC <- getWaitingTimesProb_toCancer3_store(popsC, k_cancer, store=TRUE)
waitingTimes_cancer_probs_popsD <- getWaitingTimesProb_toCancer3_store(popsD, k_cancer, store=TRUE)

#------------------------------
### waiting time to mutator ###
#------------------------------

tmp_mutatorPopsA_fn <- file.path(rdataDir,"popsA_mutatorTimes.rda")
tmp_mutatorPopsB_fn <- file.path(rdataDir,"popsB_mutatorTimes.rda")
tmp_mutatorPopsC_fn <- file.path(rdataDir,"popsC_mutatorTimes.rda")
tmp_mutatorPopsD_fn <- file.path(rdataDir,"popsD_mutatorTimes.rda")

if (!file.exists(tmp_mutatorPopsA_fn)) {
  waitingTimes_mutator_popsA <- getWaitingTimes_toMutator(popsA)
  save(waitingTimes_mutator_popsA, file=tmp_mutatorPopsA_fn)
} else load(tmp_mutatorPopsA_fn)

if (!file.exists(tmp_mutatorPopsB_fn)) {
  waitingTimes_mutator_popsB <- getWaitingTimes_toMutator(popsB)
  save(waitingTimes_mutator_popsB, file=tmp_mutatorPopsB_fn)
} else load(tmp_mutatorPopsB_fn)

if (!file.exists(tmp_mutatorPopsC_fn)) {
  waitingTimes_mutator_popsC <- getWaitingTimes_toMutator(popsC)
  save(waitingTimes_mutator_popsC, file=tmp_mutatorPopsC_fn)
} else load(tmp_mutatorPopsC_fn)

if (!file.exists(tmp_mutatorPopsD_fn)) {
  waitingTimes_mutator_popsD <- getWaitingTimes_toMutator(popsD)
  save(waitingTimes_mutator_popsD, file=tmp_mutatorPopsD_fn)
} else load(tmp_mutatorPopsD_fn)

# waitingTimes_mutator_popsB <- getWaitingTimes_toMutator(popsB)
# waitingTimes_mutator_popsC <- getWaitingTimes_toMutator(popsC)
# waitingTimes_mutator_popsD <- getWaitingTimes_toMutator(popsD)

##############
### Diversity
##############

mut4_weighted_popsA <- get4MutationClasses3_store(popsA, weighted=TRUE, store=TRUE)
mut4_weighted_popsB <- get4MutationClasses3_store(popsB, weighted=TRUE, store=TRUE)
mut4_weighted_popsC <- get4MutationClasses3_store(popsC, weighted=TRUE, store=TRUE)
mut4_weighted_popsD <- get4MutationClasses3_store(popsD, weighted=TRUE, store=TRUE)

###########################################################
###  ####

### blablabla ####

### blabla

## bla
