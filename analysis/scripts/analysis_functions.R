###############################
### pre-processing functions

require(compiler)

#' Read in a population file.
#' @param tab population file (usuallly of form r*.pop), tab separated
#' table sorted according to count of genotype
#' @param version integer element (0,1,2) specifying which implementation
#' of reading in tables should be used. 0 is standard, but optimized read.csv,
#' 1 is using a sql stage, and 2 is using data.table via fread.
readPop <- function(tabFile, version=1) {
  if (version == 0) {
    tab <- read.csv(tabFile, sep="\t", colClasses=c(
      c("integer",
        "numeric",
        "integer",
        "numeric",
        "character",
        "integer",
        "integer",
        "integer",
        "integer"))) 
    
  } else if (version == 1) {
    require(sqldf)
    f <- file(tabFile)
    tab <- sqldf("select * from f", dbname = tempfile(), 
                 file.format = list(header = T, row.names = F, sep="\t"),
                 method=c("integer",
                          "numeric",
                          "integer",
                          "numeric",
                          "character",
                          "integer",
                          "integer",
                          "integer",
                          "integer")
    )
  } else if (version == 2) {
    require(data.table)
    tab <- fread(tabFile, sep="\t")
    tab <- as.data.frame.matrix(tab)
  } else stop("version needs to be element of 0:2")
  res <- tab[order(tab[,"count"], decreasing=T),]
  res
}

readPop <- cmpfun(readPop)

switchPops2 <- function(llt) {
  if (! length(names(llt)) > 0) {
    warning("length of names(llt) <= 0")
    return(NULL)
  }
  
  firstChar <- substr(names(llt)[1], start=1, stop=1)
  
  if (length(firstChar) == 0 || is.na(firstChar))
    return(NULL)
  
  if (firstChar == 'g') {
    genNames <- names(llt)
    if (length(names(llt)) <= 0) {
      warning("length of names(llt) <= 0")
      return(NULL)
    }
#     print(genNames)
    repNames <- names(llt[[1]])
#     print(repNames)
    tmp1 <- lapply(repNames, function(rep) {
      tmp2 <- lapply(genNames, function(gen) llt[[gen]][[rep]])
      names(tmp2) <- genNames
      tmp2
    })
    names(tmp1) <- repNames
    return(tmp1)
  } else if (firstChar == 't') {
    repNames <- names(llt)
#     print(repNames)
    if (length(names(llt)) <= 0) {
      return(NULL)
      warning("length of names(llt) <= 0")
    }
    genNames <- names(llt[[1]])
#     print(genNames)
    tmp1 <- lapply(genNames, function(gen) {
      tmp2 <- lapply(repNames, function(rep) llt[[rep]][[gen]])
      names(tmp2) <- repNames
      tmp2
    }) 
    names(tmp1) <- genNames
    tmp1
  } else stop("names(llt) not properly named")
}

switchPops2 <- cmpfun(switchPops2)
switchPops3 <- function(lllt) lapply(lllt, switchPops2)

#' read replicates of a simulation stored in one dir
#' @param dir directory containing population files of format r*.pop
#' @param whatPop wildcard expression to determine which population files to be read in (e.g., "r*.pop")
#' @return tables named list of tables
readPopOfDir <- function(dir, whatPop, parallel=F, cores=3) {
  #   cwd <- getwd()
  #   setwd(dir)
  popFiles <- grep(whatPop, list.files(dir), value=T)
  popFiles <- popFiles
  if (!parallel)
    tables <- lapply(sapply(popFiles, function(x) file.path(dir,x)), readPop, version=0)
  else {
    require(parallel)
    if (length(popFiles) > 0)
      tables <- mclapply(sapply(popFiles, function(x) file.path(dir,x)), readPop, version=0,
                         mc.cores=cores, mc.preschedule=FALSE)
    else {
      warning(sprintf("dir %s empty", dir))
      return(NULL)
    }
  }
  names(tables) <- sub("(.+)\\.pop","\\1", popFiles) #sub("[a-zA-Z0]*([1-9]+)[a-zA-Z]*","\\1", popFiles)
  #   setwd(cwd)
  tables
}

#' read replicates of several simulations stored separate folders
#' @param dirs vector of directores each of which contains pop files
#' @return nested named list of named lists of population tables
readPopOfDirs <- function(dirs, whatPop, parallel=F, cores=3) {
  if (!parallel)
    tablesz <- lapply(dirs, readPopOfDir, whatPop)
  else {
    require(parallel)
    tablesz <- mclapply(dirs, readPopOfDir, whatPop, mc.cores=cores, mc.preschedule=FALSE)
  }
  names(tablesz) <- sub(".+[/\\](.+)$", "\\1", dirs)
  tablesz
}

#' partition a table, get either upper, middle, or lower part
#'♡@param tab table to be partitioned
#' @param whichPart which part of the table to be retained
#' @param size number of rows of retained part of table
getPart <- function(tab, whichPart, size=100) {
  if (whichPart == "upper") {
    i <- 1
    j <- i+size-1
  } else if (whichPart == "middle") {
    i <- (nrow(tab)/2) - floor(size/2) + 1
    j <- (nrow(tab)/2) + ceiling(size/2)
  } else if (whichPart == "lower") {
    i <- nrow(tab) - size + 1
    j <- nrow(tab)
  } else stop("'whichPart' should be either 'upper', 'middle', or ' lower'")
  
  tab[i:j,]
}

getParts1 <- function(tabs, whichPart, size=100) 
  lapply(tabs, function(tab) getPart(tab, whichPart, size))

getParts2 <- function(listOfListsOfTabs, whichPart, size=100) 
  lapply(listOfListsOfTabs, function(tabs) getParts1(tabs, whichPart, size))

getParts3 <- function(listOfListsListsOfTabs, whichPart, size=100) 
  lapply(listOfListsListsOfTabs, function(x) getParts2(x, whichPart, size))

#' @param obj list of replicates (final, e.g., "r011.pop", "r021.pop")
truncateObjectFinal <- function(obj, thresh) {
  if (length(obj) >= thresh)
    return(obj[1:thres])
  else {
    warning("cannot truncate, object too small")
    return(obj)
  }
} 

#' @param obj generation list of replicate lists (intermediate, e.g., list("g0001"=list("r010.pop")))
truncateObjectInter <- function(obj, thresh, parallel=goParallel, cores=nCores) {
  if (parallel) {
    require(parallel)
    mclapply(obj, "[", 1:thresh, mc.cores=nCores)
  } 
  else
    return(lapply(obj, "[", 1:thresh))
} 

truncateObjectInter2 <- function(lobj, thresh, parallel=goParallel, cores=nCores) {
  if (parallel) {
    require(parallel)
    return(mclapply(lobj, truncateObjectInter, thresh, mc.cores=cores)) 
  }
  else
    return(lapply(lobj, truncateObjectInter, thresh))
}

#' give a list of tables of intermediate populations of n replicates and g generations, separate
#' tables according to each generation (option 'g'), or each replicate (option 'r')
#' @param tabs list of (intermediate) populations. names should follow pattern "t[0-9]+[-_]g[0-9]+\.pop"
#' @param gr character either 'g' or 'r', indicating if tables should be clustered together according
#'        to generations or replicates
#' @return named list A of named lists B of named tables C, where each element of A contains a list
#'        of either each generation, or each replicate (depending on gr parameter)
separateIntermediatePops <- function(tabs, gr) {
  if (!all( sapply(names(tabs), function(x) substring(x,1,1)) == "t" ))
    warning("not all table names start with 't', probably not pure intermediate populations or empty?")
  
  if (gr == 'g') 
    grs_unique <- unique(sub("^t.+[-_](g.+)$", "\\1", names(tabs))) # all existing generations
  else if (gr == 'r')
    grs_unique <- unique(sub("^(t.+)[-_]g.+$", "\\1", names(tabs))) # all existing replicates
  else stop ("gr parameter needs to be either 'g' or 'r")
  
  res <- lapply(grs_unique, function(grlcl) {
    tmpNames <- grep(grlcl, names(tabs))
    lclRes <- tabs[tmpNames]
    if (gr == 'g')
      names(lclRes) <- sub("^(t.+)[-_]g.+$", "\\1", names(lclRes))
    else
      names(lclRes) <- sub("^t.+[-_](g.+)$", "\\1", names(lclRes))
    lclRes
  })
  names(res) <- grs_unique
  res
}

separateIntermediatePops <- cmpfun(separateIntermediatePops)

separateIntermediatePops2 <- function(listOfListsOfTabs, gr) 
  lapply(listOfListsOfTabs, separateIntermediatePops, gr)

#########################
### analysis functions

#' get stat values of specific columns for a set of tables
#' @param tabs list of tables containing same number of rows (nrow(table) == n)
#' @param whichColumns vector of column names contained in every table, size m
#' @return single table of size n x m containing mean values for same rows across all tables
getStatsOf <- function(tabs, f,
                       whichColumns=c(
                         "count",
                         "mean.count",
                         "number.mutations",
                         "fitness",
                         "number.drivers",
                         "number.superdrivers",
                         "number.mutators",
                         "number.passengers")) 
{
  if (length(tabs) < 1) {
    warning("length of tabs is < 1")
    return(NULL)
  } 
  if (!all(nrow(tabs[[1]]) == sapply(tabs, nrow)))
    stop("all tables are required to have equal number of rows")
  
  sapply(whichColumns, function(column) {
    tmp_table <- sapply(tabs, function(tab) {
      #       print(colnames(tab))
      tab[,column]
    })
    f(tmp_table)
  })
}

getStatsOf <- cmpfun(getStatsOf)

getMeansOf <- function(tab, ...) getStatsOf(list(tab), rowMeans, ...)

getMeansOf1 <- function(tabs, ...) getStatsOf(tabs, rowMeans, ...)

#' perform getMeansOf for a list of lists of tables
getMeansOf2 <- function(listOfListsOfTabs, ...)
  lapply(listOfListsOfTabs, getMeansOf1, ...)

#' perform getMeansOf for a list of lists of tables
getMeansOf3 <- function(listOfListsOfListsOfTabs, ...)
  lapply(listOfListsOfListsOfTabs, getMeansOf2, ...)

getSDsOf1 <- function(tab, ...) getStatsOf(tab, function(x) apply(x, 1, sd), ...)

getSDsOf2 <- function(listOfListsOfTabs, ...) 
  lapply(listOfListsOfTabs, getSDsOf1, ...)

getSDsOf3 <- function(listOfListsOfListsOfTabs, ...) 
  lapply(listOfListsOfListsOfTabs, getSDsOf2, ...)

#' perform weighted average for every column of table, where weit is relative count
#' @param tab table
#' @param col name of column that shoul be used as weight
#' @return named vector of length length(tab[,sapply(tab,is.numeric)]) of weighted averages
getWeightedAverage <- function(tab, col="mean.count") {
  we <- tab[,col]
  if (any(we < 0))
    stop("weights contain negative values!")
  numericCols <- sapply(tab, is.numeric) #apply(tab[,-5], 2, function(x) print(str(x)))
  sumwe <- sum(we)
  apply(tab[,numericCols], 2, function(x) sum(x*we)/sumwe)
}

getWeightedAverage <- cmpfun(getWeightedAverage)

getWeightedAverageOf1 <- function(tabs, col="mean.count") {
  lapply(tabs, getWeightedAverage, col)
}

getWeightedAverageOf2 <- function(listOfListsOfTabs, col="mean.count") {
  lapply(listOfListsOfTabs, getWeightedAverageOf1, col)
}

getWeightedAverageOf3 <- function(listOflistsOfListsOfTabs, col="mean.count") {
  lapply(listOflistsOfListsOfTabs, getWeightedAverageOf2, col)
}

getWeightedVariance <- function(tab, col="mean.count") {
  we <- tab[,col]  
  if (any(we < 0))
    stop("weights contain negative values!")
  numericCols <- sapply(tab, is.numeric)
  apply(tab[,numericCols], 2, function(x) {
    mu <- sum(x*we)/sum(we)
    sum(we*(x-mu))/sum(we)
  })
}

getWeightedVariance <- cmpfun(getWeightedVariance)

getWeightedVarianceOf1 <- function(tabs, col="mean.count") {
  lapply(tabs, getWeightedVariance, col)
}

getWeightedVarianceOf2 <- function(listOfListsOfTabs, col="mean.count") {
  lapply(listOfListsOfTabs, getWeightedVarianceOf1, col)
}

getWeightedVarianceOf3 <- function(listOflistsOfListsOfTabs, col="mean.count") {
  lapply(listOflistsOfListsOfTabs, getWeightedVarianceOf2, col)
}

### get one clone ###

#' Get a single clone entry from a population of clones
#' @param tab data.frame a population table
#' @param which character string specifies what type of clone should be computed
get1Clone <- function(tab, which, ...) {
  if (which == "wav") {
    return(getWeightedAverage(tab, ...))
  } else if (which == "wsd") {
    return(getWeightedVariance(tab, ...))
  } else if (which == "highestFreq") {
    return(tab[1,])
  } else if (which == "lowestFreq") {
    return(tab[nrow(tab),])
  } else stop("'which' argument needs to be element of c(\"wav\", \"wsd\", \"highestFreq\", \"lowestFreq\")")
}

get1Clone1 <- function(tabs, which, ...) lapply(tabs, get1Clone, which, ...)
get1Clone2 <- function(tabss, which, ...) lapply(tabss, get1Clone1, which, ...)
get1Clone3 <- function(tabsss, which, ...) lapply(tabsss, get1Clone2, which, ...)

### diversity ###

#' create an 4-d array of class "4mut", which indicates how many genotypes have i j k l mutations
#' @param tab population table
#' @param weighted logical if counts should be weighted. weights are not normalized, normalize before!
#' @param arrSize size of each dimension
#' @param cols chr vecor indicating which columns of tab should be used for mutation count
#' @return 4-d array, where each position i j k l returns number of genotypes with
#' i-1 j-1 k-1 and l-1 mutations
get4MutationClasses <- function(tab, weighted, arrSizes=c(100,10,10,100),
                                cols=c("number.drivers",
                                       "number.superdrivers",
                                       "number.mutators",
                                       "number.passengers"),
                                whichWeight="mean.count") {
  #   stab <- tab[,cols]
  #   cols <- c(whichWeight,cols)
  #   dimnames <- lapply(arrSizes, function(x) as.character(1:x))
  #   names(dimnames) <- cols
  #   print((dimnames))
  arr <- array(0, dim= arrSizes+1, dimnames=dimnames) #rep(arrSize, length(cols)))
  #   print(dimnames(arr))
  we <- tab[,whichWeight]
  sumwe <- sum(we)
  for (r in 1:nrow(tab)) {
    tmpVec <- unlist(tab[r,cols]) + 1 # because cannot be indexed by starting 0
    i <- tmpVec["number.drivers"] # drivers
    j <- tmpVec["number.superdrivers"] # superdrivers
    k <- tmpVec["number.mutators"] # mutators
    l <- tmpVec["number.passengers"] # passengers
    add <- ifelse(weighted, 1*we[r], 1)
    arr[i,j,k,l] <- arr[i,j,k,l] + add
  }
  class(arr) <- "4mut"
  arr
}

get4MutationClasses <- cmpfun(get4MutationClasses)

get4MutationClasses1 <- function(tabs, weighted, arrSizes=c(100,10,10,100),
                                 cols=c("number.drivers",
                                        "number.superdrivers",
                                        "number.mutators",
                                        "number.passengers"),
                                 whichWeight="mean.count",
                                 deep_parallel=F,
                                 deep_cores=1) {
  print(sprintf("We are in get4MutationClasses1 - our whichWeight is %s and the number of tabs is %s", whichWeight, length(tabs)))
  if (!deep_parallel) {
    print("We are NOT doing deep_parallel")
    res <- lapply(tabs, get4MutationClasses, weighted, arrSizes, cols, whichWeight)
  } else {
    require(parallel)
    print(sprintf("We are doing deep_parallel with %s cores", deep_cores))
    res <- mclapply(tabs, get4MutationClasses, weighted, arrSizes, cols, whichWeight, mc.cores=deep_cores)
    print("Done with get4MutationClasses1")
  }
  res
}

get4MutationClasses_mean1 <- function(tabs, weighted, arrSizes=c(100,10,10,100),
                                      cols=c("number.drivers",
                                             "number.superdrivers",
                                             "number.mutators",
                                             "number.passengers"),
                                      whichWeight="mean.count") {
  classes <- get4MutationClasses1(tabs, weighted, arrSizes=arrSizes, cols=cols, whichWeight=whichWeight)
  #   getMutatedSites <- function(listOfArr) {
  #     # create table where each row indicates indices at which position an array in listOfArr is mutated
  #     tmpRes <- do.call(rbind, lapply(listOfArr, function(x) which(x > 0, arr.ind=T)))
  #     unique(tmpRes)
  #   }
  #   resArr <- array(0, dim= arrSizes)
  #   apply
  res <- Reduce('+', classes) / length(classes)
  class(res) <- "4mut"
  res
}

get4MutationClasses_mean2 <- function(tabss, weighted, arrSizes=c(100,10,10,100),
                                      cols=c("number.drivers",
                                             "number.superdrivers",
                                             "number.mutators",
                                             "number.passengers"),
                                      whichWeight="mean.count", goParallel=goParallel,
                                      cores=1) {
  print(sprintf("We are in get4MutationClasses_mean2 and the number of tabss is %s", length(tabss)))
  if (!goParallel)
    lapply(tabss, get4MutationClasses_mean1, weighted, arrSizes=arrSizes, cols=cols, whichWeight=whichWeight)
  else {
    require(parallel)
    print(sprintf("Cores in get4MutationClasses_mean2: %s for %s tabss", cores, length(tabss)))
    mclapply(tabss, get4MutationClasses_mean1, weighted, arrSizes=arrSizes, cols=cols, whichWeight=whichWeight,
             mc.cores=cores)
  }
  
}

get4MutationClasses2 <- function(tabss, weighted, arrSizes=c(100,10,10,100),
                                 cols=c("number.drivers",
                                        "number.superdrivers",
                                        "number.mutators",
                                        "number.passengers"),
                                 whichWeight="mean.count",
                                 deep_parallel=F,
                                 deep_cores=1) {
    res <- lapply(tabss, get4MutationClasses1, weighted, arrSizes, cols, whichWeight, 
                  deep_parallel=deep_parallel, deep_cores=deep_cores)
    print("### ---- DONE WITH get4MutationClasses2 ---- ### ")
    res
}

get4MutationClasses3 <- function(tabsss, weighted, arrSizes=c(100,10,10,100),
                                 cols=c("number.drivers",
                                        "number.superdrivers",
                                        "number.mutators",
                                        "number.passengers"),
                                 whichWeight="mean.count") {
  lapply(tabsss, get4MutationClasses2, weighted, arrSizes, cols, whichWeight)
}

get4MutationClasses3_store <- function(lllt, weighted, arrSizes=c(100,10,10,100),
                                       cols=c("number.drivers",
                                              "number.superdrivers",
                                              "number.mutators",
                                              "number.passengers"),
                                       whichWeight="mean.count",
                                       store, parallel=goParallel, cores=nCores)
  getFancy(lllt, get4MutationClasses2, weighted=weighted, arrSizes=arrSizes,
           whichWeight=whichWeight, store=store, funName="4mutClasses", parallel=parallel, cores=cores)

# get4MutationClasses2

#' wrapper function to query "4mut" arrays
#' @param arr 4-d array of class "4mut"
#' @param vec named integer vector of length 4. elements are used for query
#' @return integer element of arr at position vec[1], etc. This element indicates
#' how many genotypes have vec["number.drivers"], vec["number.superdrivers"],...,
#' number of mutations
get4MutClassEntry <- function(arr, tab) {
  if (class(arr) != "4mut")
    stop("class(arr) needs to equal \"4mut\"")
  
  if (is.null(dim(tab))) {
    vec <- tab
    if (length(vec) != 4)
      stop("length(vec) needs to equal 4")
    
    vec <- vec+1
    i <- vec["number.drivers"]
    j <- vec["number.superdrivers"]
    k <- vec["number.mutators"]
    l <- vec["number.passengers"]
    res <- arr[i,j,k,l]
  } else {
    if (class(arr) != "4mut_tab")
      stop("class(tab) needs to equal \"4mut_tab\", if tab is table")
    res <- arr[tab+1]
  }
  res
}

get4MutClassEntries <- function(arr, thelist) {
  if (class(arr) != "4mut")
    stop("class(arr) needs to equal \"4mut\"")
  
  if (!setequal(names(thelist), c("number.drivers","number.superdrivers","number.mutators","number.passengers")))
      stop("names(thelist) should be number.drivers, etc")
      
  thelist <- lapply(thelist, function(x) x+1)
  
  drivers <- thelist[["number.drivers"]]
  superdrivers <- thelist[["number.superdrivers"]]
  mutators <- thelist[["number.mutators"]]
  passengers <- thelist[["number.passengers"]]
  arr[drivers,superdrivers,mutators,passengers]
}

getSumOfClasses <- function(arr, thelist=list(drivers=NA, superdrivers=NA, 
                                              mutators=NA, passengers=NA)) {
  if (class(arr) != "4mut")
    stop("class(arr) needs to equal \"4mut\"")
  
  dims <- dim(arr)
  
  entry_driver <- ifelse(is.na(thelist$drivers), 0:dims[1], thelist$drivers)
  entry_superdriver <- ifelse(is.na(thelist$superdrivers), 0:dims[2], thelist$superdrivers)
  entry_mutator <- ifelse(is.na(thelist$mutators), 0:dims[3], thelist$mutators)
  entry_passenger <- ifelse(is.na(thelist$passengers), 0:dims[4], thelist$passengers)

  sum(get4MutClassEntries(arr, thelist=list(number.drivers=entry_driver,
                                            number.superdrivers=entry_superdriver,
                                            number.mutators=entry_mutator,
                                            number.passengers=entry_passenger)))
}

get4MutClass_mutations <- function(arr) {
  if (class(arr) != "4mut")
    stop("class(arr) needs to equal \"4mut\"")
  
  res <- which(arr > 0, arr.ind=T) - 1
  colnames(res) <- c("number.drivers",
                     "number.superdrivers",
                     "number.mutators",
                     "number.passengers")
  class(res) <- "4mut_tab"
  res
}

#' @param vec named list with names(vec) == c("number.drivers", <etc>),
#' where elements are vectors 
make4mut_tab <- function(vec) {
  
}

#' get number of different mutation classes in "4mut" class object 
#' @param arr array object of class "4mut"
getNumber4mutClasses <- function(arr) {
  length(which(arr > 0))
}

getNumber4mutClasses1 <- function(arrs, goParallel=F, cores=1) {
  if (!goParallel)
    res <- lapply(arrs, getNumber4mutClasses)
  else {
    require(parallel)
#     print(sprintf("Number of cores in getNumber4mutClasses1: %s for %s arrs", cores, length(arrs)))
    res <- mclapply(arrs, getNumber4mutClasses, mc.cores=cores)
  }
  unlist(res)
}

getNumber4mutClasses2 <- function(arrss, ...) {
  lapply(arrss, getNumber4mutClasses1, ...)
}

getNumber4mutClasses3 <- function(arrsss) {
  lapply(arrsss, getNumber4mutClasses2)
}

##########################################

getOneValueStatForeach <- function(tab, f,
                                   whichColumns=c(
                                     "count",
                                     "mean.count",
                                     "number.mutations",
                                     "fitness",
                                     "number.drivers",
                                     "number.superdrivers",
                                     "number.mutators",
                                     "number.passengers")) {
  tab <- tab[,whichColumns]
  if (is.null(dim(tab)))
    return(NULL)
  apply(tab, 2, f)
}

getOneValueStatForeach <- cmpfun(getOneValueStatForeach)

getOneValueStatForeach_max <- function(tab,
                                       whichColumns=c(
                                         "count",
                                         "mean.count",
                                         "number.mutations",
                                         "fitness",
                                         "number.drivers",
                                         "number.superdrivers",
                                         "number.mutators",
                                         "number.passengers")) {
  getOneValueStatForeach(tab, max, whichColumns)
}

getOneValueStatForeach_max1 <- function(tabs, ...) lapply(tabs, getOneValueStatForeach_max, ...)
getOneValueStatForeach_max2 <- function(tabss, ...) lapply(tabss, getOneValueStatForeach_max1, ...)
getOneValueStatForeach_max3 <- function(tabsss, ...) lapply(tabsss, getOneValueStatForeach_max2, ...)

getOneValueStatForeach_min <- function(tab, 
                                       whichColumns=c(
                                         "count",
                                         "mean.count",
                                         "number.mutations",
                                         "fitness",
                                         "number.drivers",
                                         "number.superdrivers",
                                         "number.mutators",
                                         "number.passengers")) {
  getOneValueStatForeach(tab, min, whichColumns)
}

getOneValueStatForeach_min1 <- function(tabs, ...) lapply(tabs, getOneValueStatForeach_min, ...)
getOneValueStatForeach_min2 <- function(tabss, ...) lapply(tabss, getOneValueStatForeach_min1, ...)
getOneValueStatForeach_min3 <- function(tabsss, ...) lapply(tabsss, getOneValueStatForeach_min2, ...)

getOneValueStatForeach_median <- function(tab, 
                                          whichColumns=c(
                                            "count",
                                            "mean.count",
                                            "number.mutations",
                                            "fitness",
                                            "number.drivers",
                                            "number.superdrivers",
                                            "number.mutators",
                                            "number.passengers")) {
  getOneValueStatForeach(tab, median, whichColumns)
}

getOneValueStatForeach_median1 <- function(tabs, ...) lapply(tabs, getOneValueStatForeach_median, ...)
getOneValueStatForeach_median2 <- function(tabss, ...) lapply(tabss, getOneValueStatForeach_median1, ...)
getOneValueStatForeach_median3 <- function(tabsss, ...) lapply(tabsss, getOneValueStatForeach_median2, ...)

getOneValueStatForeach_Ngenotypes <- function(tab, 
                                              whichColumns=c(
                                                "count",
                                                "mean.count",
                                                "number.mutations",
                                                "fitness",
                                                "number.drivers",
                                                "number.superdrivers",
                                                "number.mutators",
                                                "number.passengers")) {
  getOneValueStatForeach(tab, length, whichColumns)
}

getOneValueStatForeach_Ngenotypes1 <- function(tabs, ...) lapply(tabs, getOneValueStatForeach_Ngenotypes, ...)
getOneValueStatForeach_Ngenotypes2 <- function(tabss, ...) lapply(tabss, getOneValueStatForeach_Ngenotypes1, ...)
getOneValueStatForeach_Ngenotypes3 <- function(tabsss, ...) lapply(tabsss, getOneValueStatForeach_Ngenotypes2, ...)

getFromOneValueStat <- function(vec, name) vec[name]
getFromOneValueStat1 <- function(vecs, name) lapply(vecs, getFromOneValueStat, name)
getFromOneValueStat2 <- function(vecss, name) lapply(vecss, getFromOneValueStat1, name)
getFromOneValueStat3 <- function(vecsss, name) lapply(vecsss, getFromOneValueStat2, name)

#### waiting time ###

#' This function gets the waiting time of a list of generations, where each elemen is a list of replicates.
#' Waiting time to k number of mutations of type m here is defined as that generation g, where among g's replicates
#' the median number of maximum numbers of m-type mutations for each replicate reaches >= k for the first time.
#' @param k numeric number of mutations to be reached
#' @param whichValue character string or vector, which indicates which type of mutation. whichValue should be the column
#' name within each population (replicate) file.
#' @param listOfListsOftabs result of separateIntermediatePops() on a single simulation run with replicates
#' @return named list of size 4
getTime_original <- function(listOfListsOfTabs, k, whichValue) { 
  
  if (length(k) > 1 && length(k) != length(whichValue))
    stop("if k is n-dimensional whichValue has to have same dimension")
  
  whichValue_merged <- do.call(paste0, as.list(whichValue))
  whichCols <- c("count",
                 "mean.count",
                 "number.mutations",
                 "fitness",
                 "number.drivers",
                 "number.superdrivers",
                 "number.mutators",
                 "number.passengers",
                 "merged")
  firstChar <- substr(names(listOfListsOfTabs)[1], start=1, stop=1)
  
  if(length(firstChar) == 0 || is.na(firstChar))
    return(NULL)
  
  ## replace original tabs with sum of whichValue entries, new column is named merged
  listOfListsOfTabs <- lapply(listOfListsOfTabs, function(x) lapply(x, function(y) {
    tmp <- y[,whichValue, drop=FALSE]
    cbind(y,merged=apply(tmp,1,sum))
  }))
  
  forGens <- function() {
    num_mutations_perGen <- lapply( # list of vectors, for each generation
      getFromOneValueStat2( # get only mutations, returns list of list of maxes
        getOneValueStatForeach_max2(listOfListsOfTabs, whichColumns=whichCols), # get maxes for each column
        "merged"),
      unlist)
    medianNum_mutations_perGen <- round(sapply(num_mutations_perGen, median))
    
    names(num_mutations_perGen) <- as.numeric(sub("g([0-9]+).*", "\\1", names(num_mutations_perGen)))
    
    gensWithLeastKmutations <- medianNum_mutations_perGen[which(medianNum_mutations_perGen >= k)]  
    tmpGens <- as.numeric(sub("g([0-9]+).*", "\\1", names(gensWithLeastKmutations)))
    names(gensWithLeastKmutations) <- tmpGens
    orderTmpGens <- order(tmpGens)
    if (length(gensWithLeastKmutations) > 0) {
      tmpRes <- gensWithLeastKmutations[orderTmpGens][1]
      gg <- names(tmpRes)
      print(gg)
      print(length(num_mutations_perGen))
      #       print(num_mutations_perGen[1:2])
      print(names(num_mutations_perGen))
      return(list(generation=gg,
                  value=tmpRes[[1]],
                  genReps=num_mutations_perGen[[gg]], # values of replicates in result generation
                  repGens=NULL))
    } else return(NULL)
  }
  
  forReps <- function() {
    num_mutations_perRep <- lapply(
      getFromOneValueStat2( # get only mutations
        getOneValueStatForeach_max2(listOfListsOfTabs, whichColumns=whichCols), # get maxes for each column
        "merged"),
      unlist)
    
    #     print(num_mutations_perRep)
    #     genNumerics <-  as.numeric(sub("g([0-9]+).*", "\\1", names(listOfListsOfTabs[[1]])))
    #     genNamesOrdered <- 
    firstGensWithLeastKmuts <- lapply(num_mutations_perRep, function(rep) {
      rep <- rep[rep >= k]
      orderGens <- order(rep, decreasing=FALSE)
      gensWithLeastKmuts <- rep[orderGens]
      genNamesNumeric <- as.numeric(sub("g([0-9]+).*", "\\1", names(gensWithLeastKmuts)))
      if (length(gensWithLeastKmuts) > 0)
        return(list(generation=genNamesNumeric[1],
                    value=gensWithLeastKmuts[[1]])) 
      else return(NULL)
    })
    
    
    firstGensWithLeastKmuts <- Filter(function(x) !is.null(x), firstGensWithLeastKmuts)
    #     firstGensWithLeastKmuts_order <- order(sapply(firstGensWithLeastKmuts, function(x) {
    # #       if (is.null(x))
    # #         return(NULL)
    #       x[["generation"]]
    #     }))
    #     firstGensWithLeastKmuts <- firstGensWithLeastKmuts[firstGensWithLeastKmuts_order]
    
    firstGensWithLeastKmuts_gensOnly <- unlist(sapply(firstGensWithLeastKmuts, function(x) {
      if (is.null(x))
        return(NULL)
      x[["generation"]]
    }))
    if (is.null(firstGensWithLeastKmuts_gensOnly))
      return(firstGensWithLeastKmuts_gensOnly)
    
    #     print(firstGensWithLeastKmuts_gensOnly)
    medianGen <- 
      as.character(sort(firstGensWithLeastKmuts_gensOnly)[[length(firstGensWithLeastKmuts_gensOnly)/2]])
    #     medianGen <- round(median(firstGensWithLeastKmuts_gensOnly))
    #     print(names(firstGensWithLeastKmuts))
    #     print(medianGen)
    #     browser()
    tmpList <- firstGensWithLeastKmuts[[medianGen]]
    c(tmpList, 
      genReps=NULL,
      list(repGens=unlist(firstGensWithLeastKmuts_gensOnly)))
  }
  
  require(compiler)
  forGens <- cmpfun(forGens)
  forReps <- cmpfun(forReps)
  
  if (firstChar == "g")
    res <- forGens()
  else if (firstChar == "t")
    res <- forReps()
  else stop("check names of list input!!")
  return(res)
}
## new (k more-dimensional)
#' This function gets the waiting time of a list of generations, where each elemen is a list of replicates.
#' Waiting time to k number of mutations of type m here is defined as that generation g, where among g's replicates
#' the median number of maximum numbers of m-type mutations for each replicate reaches >= k for the first time.
#' @param k numeric number of mutations to be reached
#' @param whichValue character string or vector, which indicates which type of mutation. whichValue should be the column
#' name within each population (replicate) file.
#' @param listOfListsOftabs result of separateIntermediatePops() on a single simulation run with replicates
#' @return named list of size 4
getTime <- function(listOfListsOfTabs, k, whichValue) { 
  
  if (length(k) > 1 && length(k) != length(whichValue))
    stop("if k is n-dimensional whichValue has to have same dimension")
  
  whichValue_merged <- do.call(paste0, as.list(whichValue))
  whichCols <- c("count",
                 "mean.count",
                 "number.mutations",
                 "fitness",
                 "number.drivers",
                 "number.superdrivers",
                 "number.mutators",
                 "number.passengers",
                 "merged")
  firstChar <- substr(names(listOfListsOfTabs)[1], start=1, stop=1)
  
  if(length(firstChar) == 0 || is.na(firstChar))
    return(NULL)
  
  if (length(k) == 1) {
    ## replace original tabs with sum of whichValue entries, new column is named merged
    listOfListsOfTabs <- lapply(listOfListsOfTabs, function(x) lapply(x, function(y) {
      tmp <- y[,whichValue, drop=FALSE]
      cbind(y,merged=apply(tmp,1,sum))
    }))
  }  
  
  #### proporsition ###
  
  if (length(k) > 1) {
    existing <- lapply(listOfListsOfTabs, function(lot) sapply(lot, function(tab) {
      tmpTab <- tab[,whichValue]
      ## which positions in tab[,colname] equals associated k value
      tmpList <- lapply(whichValue, function(colname) which(tab[,colname] == k[which(colname==whichValue)]))
      theintersect <- Reduce(intersect, tmpList)
      ifelse(length(theintersect) > 0, TRUE, FALSE) # if any rows has combination of k
#       tmpTab_rowLogicals <- apply(tmpTab, 1, 
#                                   function(row)
#                                     row[whichValue[1]] == k[1] &
#                                     row[whichValue[2]] == k[2])
#       any(tmpTab_rowLogicals)
    }))
  }
  
  forGens <- function() {
    if (length(k) == 1) {
      num_mutations_perGen <- lapply( # list of vectors, for each generation
        getFromOneValueStat2( # get only mutations, returns list of list of maxes
          getOneValueStatForeach_max2(listOfListsOfTabs, whichColumns=whichCols), # get maxes for each column
          "merged"),
        unlist)
      medianNum_mutations_perGen <- round(sapply(num_mutations_perGen, median))
      
      names(num_mutations_perGen) <- as.numeric(sub("g([0-9]+).*", "\\1", names(num_mutations_perGen)))
      
      gensWithLeastKmutations <- medianNum_mutations_perGen[which(medianNum_mutations_perGen >= k)]  
      tmpGens <- as.numeric(sub("g([0-9]+).*", "\\1", names(gensWithLeastKmutations)))
      names(gensWithLeastKmutations) <- tmpGens
      orderTmpGens <- order(tmpGens)
      if (length(gensWithLeastKmutations) > 0) {
        tmpRes <- gensWithLeastKmutations[orderTmpGens][1]
        gg <- names(tmpRes)
        print(gg)
        print(length(num_mutations_perGen))
        #       print(num_mutations_perGen[1:2])
        print(names(num_mutations_perGen))
        return(list(generation=gg,
                    value=tmpRes[[1]],
                    genReps=num_mutations_perGen[[gg]], # values of replicates in result generation
                    repGens=NULL))
      } else return(NULL) 
    } else {
      # return true if any replicate has (kr,ks) mutations
      trueGens <- subset(existing,sapply(existing, function(gen) any(gen)))
      # compute mean among replicates
      trueGens_means <- sapply(existing, function(gen) mean(gen))
      # excluded generations in which median of replicates is below 0.5
      trueGens_means_truncated <- trueGens_means[sort(which(trueGens_means > 0.5))]
      # the first generation is the desired generation
      res_generation <- as.numeric(sub("g(.+)","\\1",names(trueGens_means_truncated)[1]))
      # its median is it value
      res_value <- as.numeric(trueGens_means_truncated[1])
      # the vector of replicates is the vector of that generation that yielded res_generation
      res_genReps <- trueGens[[names(trueGens_means_truncated)[1]]]
      return(list(generation=res_generation,
                  value=res_value,
                  genReps=res_genReps,
                  repGens=NULL))
    }
  }
  ##
  forReps <- function() {
    
    if (length(k) == 1) {
      num_mutations_perRep <- lapply(
        getFromOneValueStat2( # get only mutations
          getOneValueStatForeach_max2(listOfListsOfTabs, whichColumns=whichCols), # get maxes for each column
          "merged"),
        unlist)
      
      #     print(num_mutations_perRep)
      #     genNumerics <-  as.numeric(sub("g([0-9]+).*", "\\1", names(listOfListsOfTabs[[1]])))
      #     genNamesOrdered <- 
      firstGensWithLeastKmuts <- lapply(num_mutations_perRep, function(rep) {
        rep <- rep[rep >= k]
        orderGens <- order(rep, decreasing=FALSE)
        gensWithLeastKmuts <- rep[orderGens]
        genNamesNumeric <- as.numeric(sub("g([0-9]+).*", "\\1", names(gensWithLeastKmuts)))
        if (length(gensWithLeastKmuts) > 0)
          return(list(generation=genNamesNumeric[1],
                      value=gensWithLeastKmuts[[1]])) 
        else return(NULL)
      })
      
      firstGensWithLeastKmuts <- Filter(function(x) !is.null(x), firstGensWithLeastKmuts)
      #     firstGensWithLeastKmuts_order <- order(sapply(firstGensWithLeastKmuts, function(x) {
      # #       if (is.null(x))
      # #         return(NULL)
      #       x[["generation"]]
      #     }))
      #     firstGensWithLeastKmuts <- firstGensWithLeastKmuts[firstGensWithLeastKmuts_order]
      
      firstGensWithLeastKmuts_gensOnly <- unlist(sapply(firstGensWithLeastKmuts, function(x) {
        if (is.null(x))
          return(NULL)
        x[["generation"]]
      }))
      if (is.null(firstGensWithLeastKmuts_gensOnly))
        return(firstGensWithLeastKmuts_gensOnly)
      
      #     print(firstGensWithLeastKmuts_gensOnly)
      medianGen <- 
        as.character(sort(firstGensWithLeastKmuts_gensOnly)[[length(firstGensWithLeastKmuts_gensOnly)/2]])
      #     medianGen <- round(median(firstGensWithLeastKmuts_gensOnly))
      #     print(names(firstGensWithLeastKmuts))
      #     print(medianGen)
      #     browser()
      tmpList <- firstGensWithLeastKmuts[[medianGen]]
      c(tmpList, 
        genReps=NULL,
        list(repGens=unlist(firstGensWithLeastKmuts_gensOnly)))
    } else {
#       trueGens <- subset(existing,sapply(existing, function(gen) any(gen)))
      # get, for each replicate, those generations that have the desired (kr,ks) mutation
      trueGens <- lapply(existing, function(rep) rep[which(rep)])
      print(trueGens)
      # the first generation for each replicate is the one desired
      firstTrueGens <- sapply(trueGens, function(rep) {
        if (length(rep) > 0)
          return(names(rep[1]))
        else
          return(NULL)
      })
      print(firstTrueGens)
      # get the generations names as numeric
      generationNumericNames <- unlist(sapply(firstTrueGens, function(rep)
        as.numeric(sub("g(.+)","\\1",rep))))
      print(generationNumericNames)
      # the median value among all generations is the "true" generation
      res_value <- median(generationNumericNames)
      print(res_value)
      # round the median to a discrete generation
      res_generation <- round(res_value)
      print(res_generation)
      # the vector of all first generations (for each replicates) is result vector
      res_repGens <- generationNumericNames
      print(res_repGens)
      return(list(generation=res_generation,
                  value=res_value,
                  genReps=NULL,
                  repGens=res_repGens))
    }    
  }
  
  require(compiler)
  forGens <- cmpfun(forGens)
  forReps <- cmpfun(forReps)
  
  if (firstChar == "g")
    res <- forGens()
  else if (firstChar == "t")
    res <- forReps()
  else stop("check names of list input!!")
  return(res)
}

getWaitingTime_toCancer <- function(llt, k) getTime(llt, k, c("number.drivers", "number.superdrivers"))
getWaitingTime_toMutator <- function(llt, M, infM=infMutators) {
  if (M==infM)
    return(NULL)
  else
    return(getTime(llt, M, "number.mutators"))
} 

getFancy <- function(lllt, f, ... , store, funName, parallel=goParallel, cores=nCores) {
  require(parallel)
  if (!parallel)
    cores <- 1
  result <- mclapply(names(lllt), function(x) {
    if(store) {
      tmpObjectName <- sprintf("%s_%s",x,funName)
      if (!file.exists(file.path(rdataDir,"test")))
        dir.create(file.path(rdataDir,"test"))
      tmpFileName <- file.path(rdataDir,"test",sprintf("%s.rda",tmpObjectName))
      if (!file.exists(tmpFileName)) {
        assign(tmpObjectName, f(lllt[[x]], ...))
        save(list=c(tmpObjectName), file=tmpFileName)
        return(get(tmpObjectName))
      } else {
        warning(sprintf("file %s already existing, loading file..",tmpFileName))
        load(tmpFileName)
        return(get(tmpObjectName))
      }
    } else {
      return(f(lllt[[x]], ...))
    } 
  }, mc.cores=cores)
  names(result) <- names(lllt)
  result
}

getWaitingTimes <- function(lllt, f, ..., store, funName, parallel=goParallel, cores=nCores)
  getFancy(lllt, f, ..., store=store, funName=funName, parallel=parallel, cores=cores)

getWaitingTimes_toCancer_store <- function(lllt, k, store, parallel=goParallel, cores=nCores)
  getWaitingTimes(lllt, getWaitingTime_toCancer, k=k, store=store, funName="t_cancer", parallel=parallel, cores=cores)

getWaitingTimes_toMutator_store <- function(lllt, M, store, parallel=goParallel, cores=nCores)
  getWaitingTimes(lllt, getWaitingTime_toCancer, M=M, store=store, funName="t_mutator", parallel=parallel, cores=cores)

getWaitingTimes_toCancer <- function(listOflistOfListsOfTabs, k)
  lapply(listOflistOfListsOfTabs, getWaitingTime_toCancer, k)

getWaitingTimes_toMutator2 <- function(listOflistOfListsOfTabs, M)
  lapply(listOflistOfListsOfTabs, getWaitingTime_toMutator, M)

#' get waiting time to mutator phenotype
#' @param lllt list of list of list of tables, as returned by separateIntermediatePops2()
#' @param M numeric threshold for when the mutator phenotype is said to be reached. In the 
#' case of M being NULL, individual tresholds are infered from names(lllt) for each element
#' of lllt. Otherwise, the same threshold is globally used for all elements (all parameter combinations)
getWaitingTimes_toMutator <- function(lllt, M=NULL, parallel=goParallel, cores=nCores) {
  if (is.null(M)) {
    require(parallel)
    ms <- as.numeric(sub(".+mutator([0-9]+).+", "\\1", names(lllt)))
    names(ms) <- names(lllt)
    
    if (!parallel)
      cores <- 1
    
    res <- mclapply(names(lllt), 
                    function(x) {
                      tmpResFilename <- sprintf("%s.rda",x)
                      if (!file.exists(tmpResFilename)) {
                        tmpRes <- getWaitingTime_toMutator(lllt[[x]], ms[[x]]) 
                        assign(x,tmpRes)
                        save(list=c(x), file=tmpResFilename)
                      } else {
                        warning(sprintf("file %s exists already, loading data..", tmpResFilename))
                        load(tmpResFilename)
                      }
                      return(get(x))
                    }, mc.cores=cores)
    names(res) <- names(lllt)
  } else {
    res <- getWaitingTimes_toMutator2(lllt, M)
  }
  res
}

### probabilites ###

#' Get probability of a population value computed by estimating relative frequency of that value 
#' among weighted averages for each replicate in a set of replicates
#' @param tabs list of data.frames of population tables
#' @param which character string specifies what type of clone should be computed
getProb1 <- function(tabs, which, k) {
  nums <- unlist(getFromOneValueStat1(getOneValueStatForeach_max1(tabs), which))
  length(which(nums >= k))/length(nums)
}

getProb2 <- function(tabss, which, k) lapply(tabss, getProb1, which, k)
getProb3 <- function(tabsss, which, k) lapply(tabsss, getProb2, which, k)

## waiting time probabilities for cancer phenotype

getWaitingTimesProb_toMutator1 <- function(tabs, M, infM=infMutators) {
  if (M==infM)
    return(0)
  else
    return(getProb1(tabs, "number.mutators", M))
}

getWaitingTimesProb_toMutator2 <- function(tabss, M) {
  sapply(tabss, getWaitingTimesProb_toMutator1, M)
}

getWaitingTimesProb_toMutator3 <- function(tabsss, M) {
  lapply(tabsss, getWaitingTimesProb_toMutator2, M)
}

## waiting time probabilities for mutator phenotype

getWaitingTimesProb_toCancer1 <- function(tabs, k) {
  getProb1(tabs, "number.mutations", k)
}

getWaitingTimesProb_toCancer2 <- function(tabss, k) {
  sapply(tabss, getWaitingTimesProb_toCancer1, k)
}

getWaitingTimesProb_toCancer3 <- function(tabsss, k) {
  lapply(tabsss, getWaitingTimesProb_toCancer2, k)
}

getWaitingTimesProb_toCancer3_store <- function(lllt, k, store, parallel=goParallel, cores=nCores)
  getFancy(lllt, getWaitingTimesProb_toCancer2, k=k, store=store, funName="t_cancer_probs", parallel=parallel, cores=cores)

### misc ###

shortSimulationName <- function(x) {
#   chartr("áéó", "aeo", simName)
#   gsub(".*mutator.+driver.+factor.+")
  options("scipen"=1, "digits"=10)
  x <- gsub("factor", "c", x)
  x <- gsub("mutator", "M", x)
  x <- gsub("driver", "s", x)
  x <- gsub("ratenormal", "u", x)
  options("scipen"=100, "digits"=10)
  x
}

getDriverSuperdriverCombsToK <- function(k, sizeDriver=100, sizeSuperdriver=10) {
  
}

##########################################


test <- function(dir) {
  return(c(dira = 3, dirb=4))
}
