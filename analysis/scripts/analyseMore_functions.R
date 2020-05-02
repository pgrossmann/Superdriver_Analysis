################################################################################################
### this file contains functions to perform analysis on pre-analysed data, such as waiting time 
### or mutation classes. most functions are meant to help in plotting results for a set of
### simulation results, where exactly one parameter is varied. wrapper functions exists that
### perform those analysis on all combinations of parameter variations.
### NOTE!! most functions assume that three the three types of mutation parameters, mutators,
### drivers, and factors, are to be analyzed
#################################################################################################

options("scipen"=100, "digits"=10)

#############################
### Waiting Time formulas ###

# my approximation, s_r is driver selection
tk_new <- function(k,s_r,N_f,N_0,u,d_r,d_s) {
  a <- log(N_f/N_0)
  b <- log((s_r*sqrt(2))/(u*(d_r+d_s)))
  c <- (log(N_f))^(3/2)
  d <- (log(N_0))^(3/2)
  
  (3*k*a*b)/(2*s_r*(c-d))
}

# Niko's approximation
tk_niko <- function(k,s_r,N_f,N_0,u,d_r,d_s) {
  a <- log(s_r/(u*(d_r+d_s)))
  b <- s_r * log(N_f*N_0)
  k * (a^2/b)
}

t_KrKs <- function(kr,ks,s_r,c,u,tk,N_f=10^9,N_0=10^6,d_r=100,d_s=10) {
  t_ks <- tk(k=ks,s_r=s_r*c,N_f=N_f,N_0=N_0,u=u,d_r=0,d_s=d_s)
  t_kr <- tk(k=kr,s_r=s_r,N_f=N_f,N_0=N_0,u=u,d_r=d_r,d_s=0)
  t_kr+t_ks
}
#############################

#' get parameter value of a simulation from from it's folder name, rda name, etc.
#' @param name string name of folder, rda name, etc.
#' @return named numeric vector of parameter values. names are "mutator", etc.
getParamsFromString <- function(name) {# as.numeric(
  sapply(c("mutator", "driver", "factor", "ratenormal"), 
         function(x) as.numeric(sub(sprintf(".*%s([0-9\\.]+).*",x), "\\1", name)))
}

getParamsFromString2 <- function(name, parametername){
  getParamsFromString(name)[parametername]
}

#' @param names string vector or list of simulation names
#' @return list of numeric vectors
getParamsFromStrings <- function(names) {
  tmp <- lapply(names, getParamsFromString)
  names(tmp) <- names
  tmp
} 

getParamsFromStrings2 <- function(names, parameternames) {
  lapply(getParamsFromStrings(names), function(x) as.numeric(x[parameternames]))
}

#' find names which have specifcied number of mutator mutations, etc...
#' @param names string vector of simulation names
#' @param mutator driver factor numeric or character specifying parameter value
#' @return character vector of matched names (names that have specified mutators, drivers, etc.)
findNames <- function(names, mutator, driver, factor, ratenormal) {
  mutatorS <- gsub("\\.", "\\\\.", mutator)
  driverS <- gsub("\\.", "\\\\.", driver)
  factorS <- gsub("\\.", "\\\\.", factor)
  rateS <- gsub("\\.", "\\\\.", ratenormal)
  regex <- sprintf(".*mutator%s.*driver%s.*factor%s.*ratenormal%s.*",mutatorS,driverS,factorS,ratenormal)
  #   regex <- gsub("\\.", "\\\\.", regex)
  grep(regex, names, value=T)
}

#' @param params named character vector, e.g. c(mutator=3, driver=1)
findNames2 <- function(names, params) {
  namesParams <- names(params)
  params <- as.character(params)
  names(params) <- namesParams
  params <- gsub("\\.", "\\\\.", params)
  print("names before grep")
  print(names)
  print("params before grep")
  print(params)
  matches <- lapply(names(params), function(x) {
    print("grepping")
    tomatch <- sprintf("%s%s",x,params[x])
    tonotmatch <- sprintf("%s.5",tomatch)
    tonotmatch2 <- sprintf("%s000",tomatch)
    matched <- names[grepl(tomatch,names)]
    matched <- matched[!grepl(tonotmatch,matched)]
    matched[!grepl(tonotmatch2,matched)]
#     grep(sprintf("%s%s",x,params[x]), names, value=T) 
  })
  print("before as.char")
  print(matches)
  print("after as.char")
  maches <- as.character(matches)
  print(matches)
  print("parameters")
  print(params)
  print("matches after reduce intersect")
  print(Reduce(intersect, matches))
  unique(Reduce(intersect, matches))
#   matches <- unlist(matches)
#   print("matches after unlist")
#   print(matches)
#   if (is.list(matches)) {
#     Reduce(intersect, matches)  
#   }    
#   else {
#     unique(matches)
#   }
}

#' get all parameter settings for a set of simulatin names
#' @param thenames character vector of simulation names
#' @return table with ncol == 3 and nrow == length(thenames)
#' of numeric parameter values- columns are named "mutator", etc.
#' Each row corresponds to one existing simulation output (its name)
getCombTable <- function(thenames) {
  Reduce(rbind, getParamsFromStrings(thenames))
}

#' vary one parameter and keep rest fixed. apply for each of this combinations an analysis function f
#' explanation: for each row of a table from getCombTable vary for each type of mutation (currently 3).
#' keep rest fixed and apply f on each row.
#' @param combTable table as from getCombTable()
#' @param f analysis function f accepting parameters objects, mutator, driver
#' @param objects list of pre-analysed objects. E.g. list of waiting time objects.
performForCombs <- function(combTable, f, ..., objects) {
  #for varying mutator, driver, factor, and normal mut. rate
  for (i in c("mutator", "driver", "factor", "ratenormal")) { 
    lclCombTable <- combTable
    lclCombTable[,i] <- ".*"
    lclCombTable <- unique(lclCombTable)
    apply(lclCombTable, 1, function(x) {
      #       tmpObj <- objects[findNames(names(objects), # get objects that has that mutators, drivers, factors
      #                                   mutator=x["mutator"],
      #                                   driver=x["driver"],
      #                                   factor=x["factor"])]
      f(objects, ...,
        mutator=x["mutator"],
        driver=x["driver"],
        factor=x["factor"],
        ratenormal=x["ratenormal"])
    })
  }
}

findNamesFromCombTable <- function(names, combTable) {
  thenames <- apply(combTable, 1, function(x) {
    pattern <- sapply(colnames(combTable), function(y) sprintf("%s%s",y,x[y])) # concat values and parameter names
    pattern <- sub("\\.","\\\\.", pattern) # escape regex special character
    gsub(" ", "", do.call(paste, as.list(pattern))) # vector to string
  })
  sapply(thenames, function(x) grep(x,names,value=T))
}


#' generic function to produce a single plot with fixed parameters, possibly plotting combinations of two
#' parameters in one plot
#' This function is like an outer function to plotList() functions. It should check input parameters.
#' @param onPlot character vector specifying names of parameters in simulation names that should be plotted in
#' variations in a single plot
#' @param onX string specifying which parameter name should be plotted on x
#' @param fix character vector specifying which parameters should be fixed
#' @param funPlotList function used to plot filtered list of simulation results
#' @param ylab string name of y-axis
#' @param funCenter function with which to be plotted point (y-value) is computed on a vector of replicate values
#' funCenter needs to also define how list of replicates is extracted from results! e.g. function(x) mean(x$repGens)
#' @param funSd function with which standard deviation is computed
#' @param plotSd logical specifying if error bars should be plotted
#' @param epsil numeric defines width of error bar
#' @param xloga log function or NULL. if !NULL then x-values are log transformed
#' @param yloga log function or NULL. ...
plotFun_varySome <- function(results, onPlot, onX, fix=c(), funPlotList, ylab, funCenter, funSd, 
                             plotSd=TRUE, epsil=0.001, xloga=NULL, yloga=NULL, plotApproximate=FALSE,
                             legendApprox="topright", approx2=FALSE, kCancer=NA, approxB=FALSE,
                             ...) {
  if (length(onPlot) > 2) 
    stop("only a maximum of 2 varying parameters per plot allowed")
  
  if (length(onX) > 1)
    stop("only a maximum of 1 varying parameter to be analyzed allowed")
  
  if ((approx2 | approxB) & is.na(kCancer))
    stop("'kCancer' cannot be NA if approx2 OR approxB are TRUE")
  
#   if (!is.null(xloga) & ! xloga%in%c("log2","log10"))
#     stop("xloga should be either NULL, or log2 or log10")
#   
#   if (!is.null(yloga) & ! yloga%in%c("log2","log10"))
#     stop("yloga should be either NULL, or log2 or log10")
  
  combTable <- getCombTable(names(results))  
  #   onX <- c( onX, setdiff( colnames(combTable), c(onX,onPlot) ) ) # ad rest of colnames(combTable) to onX
  
  print("onX")
  print(onX)
  print("onPlot")
  print(onPlot)
  print("combTable")
  print(combTable)
  print("names(results)")
  print(names(results))
  print("results")
  print(results)
  
  if (length(onX) + length(onPlot) > ncol(combTable)) {
    print("onX")
    print(onX)
    print("onPlot")
    print(onPlot)
    print("combTable")
    print(combTable)
    stop("number of parameters too high")
  }
  
  # no input parameter vectors should overlap
  if (length(intersect(onX, onPlot)) > 0)
    stop("combination of onX an onPlot not reasonable (intersect not 0)")  
    
  # to check if parameter names of results are same as inputted viao onX, etc.
  if (! setequal(colnames(combTable), c(onX,onPlot,names(fix))))
    stop("concat of onX, onPlot, and fix need to be parameter names of 'results'")
  
  if (plotApproximate && length(onPlot) > 1)
    stop("if approximated points should be plotted, onPlot should have length 1 only 
         in order to avoid overplotting")
  
#   for (j in onX) { # one plot for each varied parameter, because each varied one is also once fixed
#     xValues <- unique(combTable[,j])
#     print(xValues)
#     for (i in xValues) { # one plot extra each for all rest parameters which are fixed then
#       lclTab <- combTable[combTable[,j] == i,] # for all entries with value i of parameter 
#       filteredNames <- findNamesFromCombTable(names(results), lclTab)
#       fixed <- c(i)
#       names(fixed) <- j
#       xParam <- setdiff(onX,j)
# 
#       fun(results[filteredNames], xParam=xParam, fixedParams=c(fix,fixed), onPlot=onPlot)
#     }
#   }
  
  fixedNames <- findNames2(names(results), fix) ## EDIT AUGUST 3, 2018. ADDED AS.CHARACTER
  
  funPlotList(results[fixedNames], xParam=onX, onPlot=onPlot, fixed=fix, ylab=ylab, 
              funCenter=funCenter, funSd=funSd, epsil=epsil, plotSd=plotSd, xloga=xloga, 
              yloga=yloga, plotApproximate=plotApproximate, legendApprox=legendApprox, 
              approx2=approx2, kCancer=kCancer, approxB=approxB, ...)
}

#' function to plot points with lines when all values are given as vectors. also can plot error bars
#' @param xxx numeric vector of x values
#' @param yyy numeric vector of y values
#' @param sss numeric vector of standard deviation errors
#' @param colcol integer or character specifying which color to be used (corresponds to one of onX)
#' @param ltylty integer or character specifying which type of line to be used (corresponds to other of onX)
#' @param epsil numeric width of error bar
plotLinesPoints <- function(xxx, yyy, sss, colcol, ltylty, epsil=0.001, plotSd=TRUE, ...) {
  idx <- order(xxx) #Edit: Dec 17, 2016
  lines(xxx[idx], yyy[idx], col=colcol, type="o", lty=ltylty, ...)
#   axis(1, at = xxx, las=2)
  print("colcol")
  print(colcol)
  print("ltylty")
  print(ltylty)
  print("xxx")
  print(xxx)
  print("yyy")
  print(yyy)
  if (plotSd) {
    print("inSd")
    epsilon = epsil
    segments(xxx, yyy-sss,xxx, yyy+sss, col=colcol, lty=ltylty)
    segments(xxx-epsilon,yyy-sss,xxx+epsilon,yyy-sss, col=colcol, lty=ltylty)
    segments(xxx-epsilon,yyy+sss,xxx+epsilon,yyy+sss, col=colcol, lty=ltylty)
  }  
}

#' used
plotList_LinesPoints <- function(results, xParam, fixedParams, onPlot, ylab, funCenter, funSd, plotSd, 
                                 epsil=0.001, plotApproximate=FALSE, legendApprox="topright", approx2=FALSE, 
                                 kCancer=NA, approxB=FALSE, ...) {
  plotList(results, xParam=xParam, fixedParams=fixedParams, onPlot=onPlot, centerFun=funCenter, ylab=ylab, 
           plotFun=plotLinesPoints, sdFun=funSd, epsil=epsil, plotSd=plotSd, plotApproximate=plotApproximate, 
           legendApprox=legendApprox, approx2=approx2, kCancer=kCancer, approxB=approxB, ...)
}

######################################################################################################
## plot_varySome on full progression here
######################################################################################################

#' used
plot_varySome.waitingTime <- function(results, onPlot, onX, fix=c(), ylab, funCenter,
                                      funSd, plotSd=TRUE, epsil=0.001, xloga=NULL, 
                                      yloga=NULL, plotApproximate=FALSE, legendApprox="topright") {
  plotFun_varySome(results, onPlot=onPlot, onX=onX, fix=fix, funPlotList=plotList_LinesPoints, 
                   ylab=ylab, funCenter=funCenter, funSd=funSd, plotSd=TRUE, epsil=epsil,
                   xloga=xloga, yloga=yloga, plotApproximate=plotApproximate, legendApprox=legendApprox)
}

#' used
plot_varySome.waitingTimeProbs <- function(results, onPlot, onX, fix=c(), ylab, epsil=0.001,
                                           xloga=NULL, yloga=NULL, plotApproximate=FALSE,
                                           legendApprox="topright") {
  plotFun_varySome(results, onPlot=onPlot, onX=onX, fix=fix, funPlotList=plotList_LinesPoints, 
                   ylab=ylab, funCenter=function(x) x, funSd=function(x) 0, plotSd=FALSE, 
                   epsil=epsil, xloga=xloga, yloga=yloga, plotApproximate=plotApproximate,
                   legendApprox=legendApprox)
}

## waiting time cancer ##

#' potentially used
plot_varySome.waitingTime.cancer_forGens <- function(results, onPlot, onX, k, fix=c(), epsil=0.001,
                                                     xloga=NULL, yloga=NULL, plotApproximate=FALSE,
                                                     legendApprox="topright") {
  
  centerFunFun <- function(x) median(x$gensRep)
  sdFunFun <- function(x) sd(x$gensRep)
  
  plotFun_varySome(results, onPlot=onPlot, onX=onX, fix=fix, funPlotList=plotList_LinesPoints, 
                   ylab=(text=sprintf("T[k==%s]",k)), funCenter=centerFunFun, funSd=sdFunFun,
                   plotSd=TRUE, epsil=epsil, xloga=xloga, yloga=yloga, plotApproximate=plotApproximate,
                   legendApprox=legendApprox)
}

#' definitely used
plot_varySome.waitingTime.cancer_forReps <- function(results, onPlot, onX, k, fix=c(), epsil=0.001, 
                                                     xloga=NULL, yloga=NULL, plotApproximate=FALSE,
                                                     legendApprox="topright", approx2=FALSE, kCancer=NA) {
  
  centerFunFun <- function(x) median(x$repGens)
  sdFunFun <- function(x) sd(x$repGens)

  plotFun_varySome(results, onPlot=onPlot, onX=onX, fix=fix, funPlotList=plotList_LinesPoints, 
#                    ylab=sprintf("%s(%s)",expression("t"["k"]),k), funCenter=centerFunFun, funSd=sdFunFun,
                   ylab=parse(text=sprintf("T[k==%s]",k)), funCenter=centerFunFun, funSd=sdFunFun,
                   plotSd=TRUE, epsil=epsil, xloga=xloga, yloga=yloga, plotApproximate=plotApproximate,
                   legendApprox=legendApprox,approx2=approx2, kCancer=kCancer)
}

#' definitely used
plot_varySome.waitingTime_2d.cancer_forReps <- function(results, onPlot, onX, fix=c(), epsil=0.001, 
                                                        xloga=NULL, yloga=NULL, plotApproximate=FALSE,
                                                        legendApprox="topright", approxB=FALSE, 
                                                        kCancer=NA, ...) {
  
  centerFunFun <- function(x) median(x$repGens)
  sdFunFun <- function(x) sd(x$repGens)
  
  plotFun_varySome(results, onPlot=onPlot, onX=onX, fix=fix, funPlotList=plotList_LinesPoints, 
                   #                    ylab=sprintf("%s(%s)",expression("t"["k"]),k), funCenter=centerFunFun, funSd=sdFunFun,
                   ylab= sprintf("Generations to k=%s, l=%s", kCancer["ks"], kCancer["kr"]), funCenter=centerFunFun,
                   funSd=sdFunFun, plotSd=F, epsil=epsil, xloga=xloga, yloga=yloga, 
                   plotApproximate=plotApproximate, legendApprox=legendApprox, approxB=approxB,
                   kCancer=kCancer, ...)
}

## waiting time mutator ##

#' potentially used
plot_varySome.waitingTime.mutator_forReps <- function(results, onPlot, onX, fix=c(), epsil=0.001,
                                                      xloga=NULL, yloga=NULL, plotApproximate=FALSE,
                                                      legendApprox="topright", plotSd=FALSE) {
  centerFunFun <- function(x) {
    tryCatch(median(x$repGens),
             error=function(e) browser())
    
    }
  sdFunFun <- function(x) sd(x$repGens)
  
  plotFun_varySome(results, onPlot=onPlot, onX=onX, fix=fix, funPlotList=plotList_LinesPoints, 
                   ylab=expression("S"["M"]), funCenter=centerFunFun, funSd=sdFunFun,
                   plotSd=plotSd, epsil=epsil, xloga=xloga, yloga=yloga, plotApproximate=plotApproximate,
                   legendApprox=legendApprox)
}

#' potentially used
plot_varySome.waitingTime.cancer_forGens <- function(results, onPlot, onX, fix=c(), epsil=0.001,
                                                     xloga=NULL, yloga=NULL, plotApproximate=FALSE,
                                                     legendApprox="topright") {
  
  centerFunFun <- function(x) median(x$gensRep)
  sdFunFun <- function(x) sd(x$gensRep)
  
  plotFun_varySome(results, onPlot=onPlot, onX=onX, fix=fix, funPlotList=plotList_LinesPoints, 
                   ylab=expression("S"["M"]), funCenter=centerFunFun, funSd=sdFunFun,
                   plotSd=TRUE, epsil=epsil, xloga=xloga, yloga=yloga, plotApproximate=plotApproximate,
                   legendApprox=legendApprox)
}

### this block actually only loops through parameters. 
### was meant to maybe modify parameters first in case of plotting either final or intermediate pops ###

#' used
plot_varySome.waitingTime.mutator <- function(results, onPlot, onX, fix=c(), ylab, funCenter,
                                              funSd, plotSd=TRUE, epsil=0.001, xloga=NULL, yloga=NULL,
                                              plotApproximate=FALSE, legendApprox="topright") {
  plotFun_varySome(results, onPlot=onPlot, onX=onX, fix=fix, funPlotList=plotList_LinesPoints, 
                   ylab=ylab, funCenter=funCenter, funSd=funSd, plotSd=FALSE, epsil=epsil,
                   xloga=xloga, yloga=yloga, plotApproximate=plotApproximate, legendApprox=legendApprox)
}

plot_varySome.N4mutClasses <- function(results, onPlot, onX, fix=c(), ylab, funCenter,
                                       funSd, plotSd=TRUE, epsil=0.001, xloga=NULL, 
                                       yloga=NULL, plotApproximate=FALSE, legendApprox="topright") {
  plotFun_varySome(results, onPlot=onPlot, onX=onX, fix=fix, funPlotList=plotList_LinesPoints, 
                   ylab=ylab, funCenter=funCenter, funSd=funSd, plotSd=plotSd, epsil=epsil,
                   xloga=NULL, yloga=NULL, plotApproximate=plotApproximate, legendApprox=legendApprox)
}

plot_varySome.OneClones <- function(results, onPlot, onX, fix=c(), ylab, funCenter,
                                    funSd, plotSd=TRUE, epsil=0.001, xloga=NULL, 
                                    yloga=NULL, plotApproximate=FALSE, legendApprox="topright") {
#   ylab <- sub("(number)\\.(.+)", "\\1 \\2", ylab)
  
  plotFun_varySome(results, onPlot=onPlot, onX=onX, fix=fix, funPlotList=plotList_LinesPoints, 
                   ylab=ylab, funCenter=funCenter, funSd=funSd, plotSd=plotSd, epsil=epsil,
                   xloga=xloga, yloga=yloga, plotApproximate=plotApproximate, legendApprox=legendApprox)
}

######################################################################################################
## finalPops here ##
######################################################################################################

#' used
plot_varySome.N4mutClasses.finalPops <- function(results, onPlot, onX, fix=c(), epsil=0.001, plotSd=TRUE,
                                                 xloga=NULL, yloga=NULL, plotApproximate=FALSE,
                                                 legendApprox="topright") {
  centerFunFun <- function(x) mean(x)
  sdFunFun <- function(x) sd(x)
  plot_varySome.N4mutClasses(results, onPlot=onPlot, onX=onX, fix=fix, ylab="no. mutation classes [n(s_r)]",
                             funCenter=centerFunFun, funSd=sdFunFun, plotSd=plotSd, epsil=epsil,
                             xloga=xloga, yloga=yloga, plotApproximate=plotApproximate,
                             legendApprox=legendApprox)
}

##

plot_varySome.oneClonesMutations.finalPops <- 
  function(results, onPlot, onX, fix=c(), whichMutation, epsil=0.001, plotSd=TRUE, 
           xloga=NULL, yloga=NULL, plotApproximate=FALSE, legendApprox="topright") {
    centerFunFun <- function(x) mean(x)
    sdFunFun <- function(x) sd(x)
    
    ylab <- sprintf("number %s",whichMutation)
    
    plot_varySome.OneClones(results, onPlot=onPlot, onX=onX, fix=fix, ylab=ylab,
                            funCenter=centerFunFun, funSd=sdFunFun, plotSd=plotSd, epsil=epsil, 
                            xloga=xloga, yloga=yloga, plotApproximate=plotApproximate,
                            legendApprox=legendApprox)
}

###########################################

plotList_test2 <- function(results, xParam, fixedParams, onPlot, ylab, funCenter, funSd, epsil, ...) {
  plotList(results, xParam, fixedParams, onPlot, funCenter, ylab, plotLinesPoints, funSd, epsil, ...)
}

plotList_test <- function(results, xParam, fixedParams, onPlot, ylab,
                          fff=function(x) mean(x$repGens), FFF=function(x) sd(x$repGens), ...) {
  plotList(results, xParam, fixedParams, onPlot, fff, ylab,
           points, FFF, ...)
}

plotFun_varySome_test <- function(results, onPlot, onX, fix=c(), ylab, fff, FFF) {
  plotFun_varySome(results, onPlot, onX, fix, plotList_test, ylab, fff, FFF)
}

plotFun_varySome_test2 <- function(results, onPlot, onX, fix=c(), ylab, funCenter, funSd, epsil) {
  plotFun_varySome(results, onPlot=onPlot, onX=onX, fix=fix, funPlotList=plotList_test2, 
                   ylab=ylab, funCenter=funCenter, funSd=funSd, plotSd=TRUE, epsil=epsil)
}

#' generic function to plot results on a filtered list of results. this list is transformed to character
#' vectors by inputted functions such as centerFun()
#' Specifications (instances) of this function are used in plot_varySome() functions.
#' @param xParam character string specifying which parameter is the one to be plotted on x-axis
#' @param fixedParams named integer vector. names are mutator, driver, etc. elements are parameter values
#' @param onPlot character vector names of parameters to be plotted at once
#' @param centerFun function to computed center from replicates. see plotFun_varySome()
#' @param ylab label of y-axis
#' @param plotFun function to use to plot points when given as vectors
#' @param sdFun function to compute error from replicates. see plotFun_varySome()
#' @param epsil numeric defines width of error bar. see plotFun_varySome()
#' @param graphPaper quick-and-dirty parameter to manually adjust graphics for paper
plotList <- function(results, xParam, fixedParams, onPlot, centerFun, ylab, plotFun, sdFun, 
                     epsil=0.001, plotSd=TRUE, plotLegend=TRUE, filterOutCond=is.null, 
                     xloga=NULL, yloga=NULL, plotApproximate=FALSE, legendApprox="topright", approx2=FALSE, 
                     kCancer=NA, approxB=FALSE, graphPaper=TRUE, plotLegend.topRight=FALSE, orderLegend=TRUE, 
                     errorCorrection.lm=FALSE, ...) {
  ## quick and dirty: if repGens empty test switchedRepGens
#   results <- sapply(results, function(x) {
#     if (!is.null(x$repGens)) {
#       return(x)
#     } else {
#       res <- x
#       res[["repGens"]] <- res$mutatorWaitingTimeSwitched$repGens
#       return(res)
#     }
#   })
  ## this is only a quick-and-dirty filter to get rid of results that did not converge! 
  ## fix this is getTime function in parts applied to k=c(kr,ks)!!
  tryCatch(
    
  results <- Filter(function(x) length(x$repGens) != 0, results)#results[5:length(results)]
  , error=function(e) browser())
  print("inplotlist")
  results <- Filter(function(x) ! filterOutCond(x), results)
  xValues <- unlist(getParamsFromStrings2(names(results), xParam))
  if (!is.null(xloga))
    xValues <- xloga(xValues)
  xli <- range(c(0,xValues))
  yValues <- sapply(results, centerFun)
  if (!is.null(yloga))
    yValues <- yloga(yValues)
  yli <- range(yValues)
  yValues_sd <- sapply(results, sdFun)
  
  print("yValues")
  print(yValues)
  print("yValues stop")
  options("scipen"=1, "digits"=10)
  tit <- do.call(
    paste, 
      lapply(names(fixedParams), function(x) sprintf("%s=%s",x,fixedParams[x]))
    )
  options("scipen"=100, "digits"=10)
  
  options("scipen"=1, "digits"=10)
  tit <- gsub("factor", "c", tit)
  tit <- gsub("mutator", "M", tit)
  tit <- gsub("driver", "s", tit)
  tit <- gsub("ratenormal", "u", tit)
  options("scipen"=100, "digits"=10)
  
  xParam <- gsub("factor", "c", xParam)
  xParam <- gsub("mutator", "M", xParam)
  xParam <- gsub("driver", "s", xParam)
  xParam <- gsub("ratenormal", "u", xParam)
  
  # get simulation names with fixed input parameters
  fixedNames <- findNames2(names(results), fixedParams)
  
  if (plotLegend & !plotLegend.topRight)
    layout(rbind(1,2), heights=c(9,1)) 
  

  if (!is.null(xloga))
    xlab <- sprintf("log(%s)",xParam)
  if (!is.null(yloga))
    ylab <- sprintf("log(%s)",ylab)
  
  # dirty to make "driver selection" out of "s"
  if (xParam == "s" & is.null(xloga))
    xlab <- "Driver selection (s)"
  
  if (any(is.infinite(c(xli,yli)))) {
    plot(1,1,type="n",main=sprintf("sorry results NULL, %s", tit), ylab=ylab, xlab=xParam)
    return(NULL)
  }

  # just comment: las=2 does not work here because of xaxt="n"
  if (!graphPaper) {
    plot(1,1,type="n",xlim=xli,ylim=yli,main=tit,xlab=xlab,ylab=ylab,xaxt="n",las=2) 
    axis(1, at = xValues, labels=TRUE, las=2, cex.axis=.77)
    axis(1, at = c(0), las=2, tck=-0.045, labels=TRUE, lwd=3, cex.axis=.7)
  } else {
    op <- par(mar = c(7,7,4,2) + 0.1) ## default is c(5,4,4,2) + 0.1
    plot(1,1,type="n",xlim=xli,ylim=yli,main=tit,xlab=xlab,ylab=ylab, xaxt="n",las=2, 
         cex.axis=1.5, ann=F, yaxs="i", xaxs="i") 
    axis(1, at=xValues, labels=F, cex.axis=1.3)
    text(xValues, 0, labels=xValues, srt=45, pos=1, xpd=T)
    title(ylab=ylab,line=5,xlab=xlab,cex.lab=1.7)
    par("usr")
    par(op)
#     axis(1, at = xValues, labels=TRUE, las=2, cex.axis=1.3) #Edit: commented out Dec 17, 2016
#     axis(1, at = c(0), las=2, tck=-0.045, labels=TRUE, lwd=3, cex.axis=.7)
  }
    
  
  # unique values (list of vectors) for each of onPlot
  onPlotUniqueVals <- lapply(onPlot, function(x) unique(getParamsFromStrings2(names(results),x)))
  onPlotUniqueVals <- lapply(onPlotUniqueVals, unlist)
  names(onPlotUniqueVals) <- onPlot  
  
  print("bla")
  print(onPlotUniqueVals)
  
  # for legend of plotting approximations
  bestFits <- list()
  bestFitsPreds <- list()
  
  # for each combination of  unique values of parameters in onPlot
  for (i in onPlotUniqueVals[[1]]) {
    
    itsBigger <- length(onPlotUniqueVals) > 1
    if (itsBigger) # if onPlot is not only of length 1
      jValues <- onPlotUniqueVals[[2]]
    else
      jValues <- 1
      
      for (j in jValues) {
        if (itsBigger)
          tmpParamVec <- c(i,j)
        else
          tmpParamVec <- c(i)
        names(tmpParamVec) <- onPlot
        #       if(length(names(results))==0)
        tmpNames <- findNames2(names(results), tmpParamVec)
        #       tmpResults <- results[tmpNames]
        colo <- which(onPlotUniqueVals[[1]] == i)
        if (itsBigger)
          ltyp <- which(onPlotUniqueVals[[2]] == j)
        else
          ltyp <- 1
        tmpPointsX <- xValues[tmpNames]
        tmpPointsY <- yValues[tmpNames]
        tmpPoints_sd <- yValues_sd[tmpNames]
        print("tmpParamVec before plotFun")
        print(tmpParamVec)
        print("tmpNames before plotFun")
        print(tmpNames)
        print("tmpPoints before plotFun")
        print(tmpPointsY)
        
        print("before plotFun")
        plotFun(tmpPointsX, tmpPointsY, tmpPoints_sd, colcol=colo, ltylty=ltyp, epsil=epsil, plotSd=plotSd, ...)
        plottedWaitingTimes.simulated <- list(xx=tmpPointsX, yy=tmpPointsY) # EDIT: AUGUST 3, 2018
        plottedWaitingTimes.theoretical <- list(xx=NA, yy=NA)
        
        if (approxB) {
          
#           xx <- seq(0.001,0.1,by=0.001) #c(0.001,0.005,0.01,0.05,0.1) #Edit commented out Dec 17, 2016
          # xx <- c(0.005, 0.01, 0.02, 0.025, 0.03, 0.04, 0.05) #Edit commented out Dec 17, 2016
          xx <- getParamsFromStrings2(names(tmpPointsX), "driver")
          yy <- t_KrKs(kr=kCancer["kr"], ks=kCancer["ks"], s_r=as.numeric(xx), c=i, u=fixedParams["ratenormal"], tk=tk_niko)
          
          yy.orig <- yy
          if (!is.null(yloga)) yy <-yloga(yy)
          ## EDIT: AUGUST 27, 2018
          if (errorCorrection.lm) {
            fn <- "waiting-time-error-correction-models.rda"
            if (file.exists(fn)) {
              print("ERROR CORRECTION ON LM BASIS")
              load(fn)
              mod <- waitingTimeErrorCorrectionModels$mods$FactorDriverKsKr
              save(mod, file="error-correction_model-used.rda")
              covar.driver <- as.numeric(getParamsFromStrings2(names(tmpPointsX), "driver"))
              covar.factor <- as.numeric(getParamsFromStrings2(names(tmpPointsX), "factor"))
              covar.ks <- as.numeric(gsub(".*waitingTimeResKr[0-9]*Ks([0-9]*).*", "\\1", names(tmpPointsX)))
              covar.kr <- as.numeric(gsub(".*waitingTimeResKr([0-9]*)Ks[0-9]*.*", "\\1", names(tmpPointsX)))
              new_df <- data.frame(Driver=covar.driver, Factor=covar.factor, Ks=covar.ks, Kr=covar.kr)
              new_df$DriverFactorRatio <- new_df$Driver/new_df$Factor
              err <- predict(mod, new_df)
              yy_corrected <- yy+err
              lines(xx, yy_corrected, col=colo, lty=5)
            } else print("NO ERROR CORRECTION MODELS AVAILABLE")
          } else {
            lines(xx, yy, col=colo, lty=5)
          }
          
#           legend("topright",legend=c("tk_niko 1e-8","tk_niko 1e-7","tk_new 1e-8","tk_new 1e-7"),
#                  col=c("black","red","black","red"), lty=5)
          plottedWaitingTimes.theoretical <- list(xx=xx, yy=yy) # EDIT: AUGUST 3, 2018
        }
        
        # EDIT: AUGUST 3, 2018
        plottedWaitingTimes <- list(simulated=plottedWaitingTimes.simulated, theoretical=plottedWaitingTimes.theoretical)
        names(plottedWaitingTimes$theoretical$xx) <- names(plottedWaitingTimes$simulated$xx)
        names(plottedWaitingTimes$theoretical$yy) <- names(plottedWaitingTimes$simulated$yy)
        
        print("after plotFun")
        
        if (length(tmpPointsX) > 2 & length(tmpPointsY) > 2) {
          fits <- getFits(tmpPointsX, tmpPointsY)# APPROACH 1
          print(fits)# APPROACH 1
          tmpList <- list(getBestFit(fits))# APPROACH 1
          
#           fits <- getFits_symbolicRegression(tmpPointsX, tmpPointsY)# APPROACH 2
#           print(fits)# APPROACH 2
#           tmpList <- list(fits)# APPROACH 2
#           preds <- predict(fits, data.frame(x=tmpPointsX, y=tmpPointsY))
#           print(preds)
#           tmpList2 <- list(preds)
#           names(tmpList2) <- colo
#           bestFitsPreds <- c(bestFitsPreds,tmpList2)
          
          names(tmpList) <- colo
          bestFits <- c(bestFits,tmpList)
        } #else bestFit <- NULL
      } 
  }
  if (graphPaper) par(op)
  
  if (plotApproximate) {
    for (fitColor in names(bestFits)) {
      thexx <- seq(xli[1], xli[2], length=50)
      theyy <- predict(bestFits[[fitColor]], data.frame(x=thexx))# APPROACH 1
#       theyy <- predict(bestFits[[fitColor]], newdata = data.frame(z=thexx)) # gives error
#       theyy <- bestFitsPreds[[fitColor]]
#       theyy <- getBestFitFunction_symbolicRegression(bestFits[[fitColor]])(thexx)# APPROACH 2
      if (length(thexx) == length(theyy))
        lines(thexx, theyy, col=fitColor, lty=5)
    }
    if (length(bestFits) > 0) {
      legendTexxt <- sapply(bestFits, function(x) deparse(formula(x))) # APPROACH 1
#       legendTexxt <- sapply(bestFits, # APPROACH 2
#                             function(x) deparse(getBestFitFunction_symbolicRegression(x))[2])
      print("bestFits")
      print(bestFits)
      print("legendtexxt")
      print(legendTexxt)
      legend(legendApprox, legend=legendTexxt, col=names(bestFits), lty=5)
    }    
  }
  
  if (approx2) {
    xx <- seq(0.001, 0.1, by=0.005) #c(0.001,0.005,0.01,0.05,0.1)
    yy <- tk_niko(k=kCancer,s_r=xx,N_f=10^9,N_0=10^6,u=1e-8,d_r=100,d_s=10)
    lines(xx,yy,col="red",type="l",xaxt="n",lty=3)
    yy <- tk_niko(k=kCancer,s_r=xx,N_f=10^9,N_0=10^6,u=1e-7,d_r=100,d_s=10)
    lines(xx,yy,col="black",lty=3)
    ##
    xx <- seq(0.001, 0.1, by=0.005) #c(0.001,0.005,0.01,0.05,0.1)
    yy <- tk_new(k=kCancer,s_r=xx,N_f=10^9,N_0=10^6,u=1e-8,d_r=100,d_s=10)
    lines(xx,yy,col="red",lty=2)
    yy <- tk_new(k=kCancer,s_r=xx,N_f=10^9,N_0=10^6,u=1e-7,d_r=100,d_s=10)
    lines(xx,yy,col="black",lty=2)
#     axis(1, at=xx)
    legend("topright",legend=c("tk_2007 1e-8","tk_2007 1e-7","tk_new 1e-8","tk_new 1e-7"),
           col=c("red","black","red","black"), lty=c(3,3,2,2))
  }  
  
  if (plotLegend) {
    isBigger <- length(onPlotUniqueVals) > 1
    len2 <- ifelse(isBigger, length(onPlotUniqueVals[[2]]), 0)
    
    #     lg_ncol <- max(length(onPlotUniqueVals[[1]]),length(onPlotUniqueVals[[2]]))
    allLength <- length(onPlotUniqueVals[[1]])+len2
    lg_ncol <- ifelse(allLength > 5, 5, allLength)
    isSecondMax <- which.max(c(length(onPlotUniqueVals[[1]]), len2))
    
    options("scipen"=1, "digits"=10)
    lg_text_cols <- sapply(onPlotUniqueVals[[1]], function(x) paste(names(onPlotUniqueVals)[1], x, sep=" = "))
    lg_text_lins <-
      if (isBigger)
        sapply(onPlotUniqueVals[[2]], function(x) paste(names(onPlotUniqueVals)[2], x, sep=" = "))
    else
      c()
    
    lg_text <- c(lg_text_cols, lg_text_lins)
    options("scipen"=100, "digits"=10)
    
    lg_cols_cols <- 1:length(onPlotUniqueVals[[1]])
    lg_cols_lins <- 
      if (isBigger)
        rep(1, len2) # black
    else
      c()
    if (isBigger)
      lg_cols <- c(lg_cols_cols, lg_cols_lins)
    else
      lg_cols <- lg_cols_cols
    
    lg_lty_cols <- rep(0, length(onPlotUniqueVals[[1]])) # no linetype
    if (isBigger)
      lg_lty_lins <- 1:len2
    else
      lg_lty_lins <- c()
    if (isBigger)
      lg_lty <- c(lg_lty_cols, lg_lty_lins)
    else
      lg_lty <- lg_lty_cols
    
    lg_lwd <- lg_lty
    
    lg_pch_cols <- rep(22, length(onPlotUniqueVals[[1]]))
    lg_pch_lty <- rep(NA, len2)
    if (isBigger)
      lg_pch <- c(lg_pch_cols, lg_pch_lty)
    else
      lg_pch <- lg_pch_cols
    
    lg_pt.bg_cols <- lg_cols_cols
    lg_pt.bg_lins <- rep(NA, len2)
    if (isBigger)
      lg_pt.bg <- c(lg_pt.bg_cols, lg_pt.bg_lins)
    else
      lg_pt.bg <- lg_pt.bg_cols
    
    lg_text <- sapply(lg_text, function(x) {
      x <- gsub("factor", "c", x)
      x <- gsub("mutator", "M", x)
      x <- gsub("driver", "s", x)
      x <- gsub("ratenormal", "u", x)
    })
    
    if (isSecondMax) {
      lg_text <- rev(lg_text)
      lg_cols <- rev(lg_cols)
      lg_lty <- rev(lg_lty)
      lg_lwd <- rev(lg_lwd)
      lg_pch <- rev(lg_pch)
      lg_pt.bg <- rev(lg_pt.bg)
    }
    if (orderLegend) {
      iix <- order(lg_text, decreasing=F)
      lg_text <- lg_text[iix]
      lg_cols <- lg_cols[iix]
      if (length(lg_lty)==length(iix))
        lg_lty <- lg_lty[iix]
      if (length(lg_lwd)==length(iix))
        lg_lwd <- lg_lwd[iix]
      if (length(lg_pch)==length(iix))
        lg_pch <- lg_pch[iix]
      if (length(lg_pt.bg)==length(iix))
        lg_pt.bg <- lg_pt.bg[iix]
    }
    if (!plotLegend.topRight) {
      parop <- par(mar=c(0, 0, 0, 0))
      plot.new()
      legend("center", legend = lg_text, 
             bty = "n", ncol=lg_ncol,
             col = lg_cols,
             lty = lg_lty, lwd = lg_lwd,
             pch = lg_pch,
             pt.bg = lg_pt.bg,
             pt.cex = 2)
      par(parop)
    } else {
      legend("topright", legend = lg_text, 
             bty = "n", ncol=1,
             col = lg_cols,
             lty = lg_lty, lwd = lg_lwd,
             pch = lg_pch,
             pt.bg = lg_pt.bg,
             pt.cex = 3, cex=1.5)
    }
  }
  return(list(bestFits=bestFits, plottedWaitingTimes=plottedWaitingTimes)) # EDIT: AUGUST 3, 2018
#   stop("bla")
#   shuffledResults_names <- lapply(onPlot[1], # e.g. c("M","ratenormal")[1]
#                             function(x) {
#                               tmpList <- lapply(onPlot[2], # e.g. 5:4 for "M"
#                                                 function(y) {
#                                                    tmpVec <- c(x,y)
#                                                    names(tmpVec) <- names(onPlot)
#                                                    findNames2(names(results), tmpVec)
#                                                  })
#                               names(tmpList) <- onPlot[2]
#                               tmpList
#                             })
#   names(shuffledResults_names) <- onPlot[1]
#   for (i in shuffledResults_names) {
#     for (j in shuffledResults_names[[i]]) {
#       print(1)
#     }
#   }
  
}

##########################
### summary by heatmap ###
##########################

# getWaitingTimeTable <- function(case, k, getValueFun, yParam="driver") {
#   # waiting time result files
#   wtFiles <- grep(sprintf("%s.*waitingTimeResK%s.*\\.rda",case,k), list.files(intermediateDir), value=T)
#   params <- getParamsFromStrings(wtFiles)
#   print(params)
#   print(wtFiles)
#   if (length(params) > 0)
#     params_unique <- sapply(names(params[[1]]),
#                             function(x) unique(
#                               as.character(sapply(params, function(y) y[x])))) # use as.char!!
#   else
#     params_unique <- c()
#   sapply(params_unique[yParam],
# }

## value ##

getHeatmapEntry_binary <- function(obj, p=0.5) {

  
  if (is.null(obj$cancerWaitingTimeProbs))
    stop("obj shouldn't be null")
  
  if (any(obj$cancerWaitingTimeProbs >= p)) 1 else 0
}

getHeatmapEntry_highest <- function(obj, p=0.95) {
  
  if (is.null(obj$cancerWaitingTimeProbs))
    stop("obj shouldn't be null")
  
  max(obj$cancerWaitingTimeProbs)
}

getHeatmapEntry_smallestGeneration <- function(obj, p=0.95) {
  
  if (is.null(obj$cancerWaitingTimeProbs))
    stop("obj shouldn't be null")
  
  genWithGreaterEqualP <- obj$cancerWaitingTimeProbs[which(obj$cancerWaitingTimeProbs >= p)]
  
  if (length(genWithGreaterEqualP) == 0)
    return(NULL)
  else {
    return(names(genWithGreaterEqualP)[1])
  }
}

getHeatmapOverlay_smallestGeneration <- function(obj, p=0.95) {
  
  if (is.null(obj$cancerWaitingTimeProbs))
    stop("obj shouldn't be null")
  
  genWithGreaterEqualP <- obj$mutatorWaitingTimeProbs[which(obj$cancerWaitingTimeProbs >= p)]
  
  if (length(genWithGreaterEqualP) == 0)
    return(NULL)
  else {
    return(names(genWithGreaterEqualP)[1])
  }
}

## text ##

getHeatmapOverlayText_binary <- function(obj, p=0.5) {
  
#   if (is.na(val))
#     return("")
#   
#   ifelse(val == 1, "*", "")
  
  if (is.null(obj$mutatorWaitingTimeProbs))
    stop("obj shouldn't be null")
  
  ifelse(any(obj$mutatorWaitingTimeProbs >= p), "*", "")
}

## get table ##

#' @param uniqueParams list of unique parameter values
getWaitingTimeTable <- 
function(case, k, getValueFun, getOverlayTextFun, yParam="driver",
         xParams=c("mutator","factor"), fixedParam=c("ratenormal"=0.00000001),
         uniqueParams=NULL, heatmap_realValue_fun=NULL, heatmap_realText_fun=NULL) {
  
  # waiting time result files
  wtFiles <- grep(sprintf("%s.*waitingTimeResK%s.*\\.rda",case,k), list.files(intermediateDir), value=T)
  wtFiles <- findNames2(wtFiles, fixedParam)
  params <- getParamsFromStrings(wtFiles)
  ##########################################################
  combTable <- getCombTable(wtFiles)  
  #   onX <- c( onX, setdiff( colnames(combTable), c(onX,onPlot) ) ) # ad rest of colnames(combTable) to onX
  if (is.null(combTable))
    return(NULL)
  if (length(yParam) + length(xParams) + length(fixedParam) > ncol(combTable)) {
    stop("number of parameters too high")
  }
  # no input parameter vectors should overlap
  if (length(intersect(yParam, intersect(xParams, names(fixedParam)))) > 0)
    stop("combination of yParm, xParams, and fixedParam not reasonable (intersect not 0)")  
  
  # to check if parameter names of results are same as inputted viao onX, etc.
  if (! setequal(colnames(combTable), c(xParams,yParam,names(fixedParam))))
    stop("concat of yParm, xParams, and fixedParax need to be parameter names of read case files")
  ############################################################
  
  if (length(params) > 0)
    params_unique <- sapply(names(params[[1]]),
                            function(x) unique(
                              as.character(sapply(params, function(y) y[x])))) # use as.char!!
  else
    params_unique <- c()
  
  #--------------------------------
  if (!is.null(uniqueParams))
    params_unique <- uniqueParams
  #--------------------------------
  
  #   yNames <- paste0(yParam, params_unique[[yParam]])
  yParamValues <- params_unique[[yParam]]
  names(yParamValues) <- paste0(yParam, params_unique[[yParam]])
  yParamValues <- lapply(yParamValues, function(x) {
    tmpV <- c(x)
    names(tmpV) <- yParam
    tmpV
  })
  
  # list of vectors
  if (length(xParams) == 1) {
    #     xParamValues <- paste0(xParams, params_unique[xParams])
    #     xParamValues <- 
    stop("length(xParams) == 1 not implemented yet")
  }
  else if (length(xParams) == 2) { # will result list of vectors 
    xParamValues <- unlist(sapply(params_unique[[xParams[1]]], 
                                  function(x) sapply(params_unique[[xParams[2]]], 
                                                     function(y) {
                                                       tmpV <- c(x,y)
                                                       names(tmpV) <- xParams
                                                       tmpV
                                                     }, simplify=FALSE), simplify=FALSE),
                           recursive=FALSE)
    names(xParamValues) <- sapply(xParamValues, function(x){
      tmp2 <- sapply(names(x), function(y) paste0(y,x[y]))
      do.call(paste0, as.list(tmp2))
    })
    xParamValues
  }
  
  else
    stop("more than 2 xParams not yet implemented")
  
  
  if (is.null(uniqueParams))
    xParamValues <- Filter(function(x) length(findNames2(wtFiles, x)) > 0, xParamValues)
  
  # heatmap values
  tableValues <- sapply(xParamValues, 
                        function(x) sapply(yParamValues,
                                           function(y) {
                                             objName <- findNames2(wtFiles, c(x,y))
                                             tryCatch(getValueFun(
                                               loadIfNot(file.path(intermediateDir,objName))),
                                                             error = function(e) {
                                                               print(e)
#                                                                browser()
                                                               NA
                                                             })
#                                              getValueFun(obj)
                                           }))
  # text to overlay ("*)
  tableOverlayText <- sapply(xParamValues, 
                             function(x) sapply(yParamValues,
                                                function(y) {
                                                  objName <- findNames2(wtFiles, c(x,y))
                                                  tryCatch(getOverlayTextFun(
                                                    loadIfNot(file.path(intermediateDir,objName))),
                                                                  error = function(e) NA)
#                                                   getOverlayTextFun(obj)
                                                }))
  ############################
  #### real values for CSV ###
  
  if (!is.null(heatmap_realValue_fun) && !is.null(heatmap_realText_fun)) {
    # heatmap values
    tableValues_realValue <- sapply(xParamValues, 
                                    function(x) sapply(yParamValues,
                                                       function(y) {
                                                         objName <- findNames2(wtFiles, c(x,y))
                                                         tryCatch(heatmap_realValue_fun(
                                                           loadIfNot(file.path(intermediateDir,objName))),
                                                                  error = function(e) {
                                                                    print(e)
                                                                    #                                                                browser()
                                                                    NA
                                                                  })
                                                         #                                              getValueFun(obj)
                                                       }))
    # text to overlay ("*)
    tableOverlayText_realValue <- sapply(xParamValues, 
                                         function(x) sapply(yParamValues,
                                                            function(y) {
                                                              objName <- findNames2(wtFiles, c(x,y))
                                                              tryCatch(heatmap_realText_fun(
                                                                loadIfNot(file.path(intermediateDir,objName))),
                                                                       error = function(e) NA)
                                                              #                                                   getOverlayTextFun(obj)
                                                            }))
  } else {
    tableValues_realValue <- NULL
    tableOverlayText_realValue <- NULL
  }
  
#   print(tableValues)
#   browser()
#   stop("asdf")
  
  #   if (length(xParams) == 1)
  #     xNames <- paste0(xParams, params_unique[xParams])
  #   else if (length(xParams) == 2) {
  #     xNames <- as.vector(unlist(sapply(params_unique[[xParams[1]]], 
  #                                       function(x) unlist(sapply(params_unique[[xParams[2]]], 
  #                                                                 function(y) paste0(xParams[1], x,
  #                                                                                    xParams[2], y))),
  #                                       simplify=FALSE)))
  #     xNames <- as.vector(unlist(sapply(params_unique[[xParams[1]]], 
  #                                       function(x) unlist(sapply(params_unique[[xParams[2]]], 
  #                                                                 function(y) {
  #                                                                   tmpV <- c(x,y)
  #                                                                   names(tmpV) <- xParams
  #                                                                   findNames2(wtFiles, tmpV)
  #                                                                 })),
  #                                       simplify=FALSE)))
  
  #     xNames <- lapply(params_unique[[xParams[1]]], 
  #                                       function(x) unlist(sapply(params_unique[[xParams[2]]], 
  #                                                                 function(y) {
  #                                                                   tmpV <- c(x,y)
  #                                                                   names(tmpV) <- xParams
  #                                                                   tmpV
  #                                                                 })))
  #                                                                 
  #     
  #   }
  #     
  #   else
  #     stop("more than 2 xParams not yet implemented")
  
  #   print(xNames)
  #   stop("asdf")
  #   sapply(xNames, function(x) print(grep(paste0(x,yNames[1]), wtFiles, value=T)))
  
  #   sapply(params_unique[yParam], function(y) sapply(xParams, function(x)
  #     do.call(paste0, )))
  
  #   yNames <- paste0()
  
  #   sapply(params_unique[yParam],)
  
  list(tableValues=tableValues,
       tableOverlayText=tableOverlayText,
       tableValues_realValue=tableOverlayText_realValue,
       tableOverlayText_realValue=tableOverlayText_realValue)
  
}

getColorCode <- function(strings, colors=c("A"="pink","B"="blue","C"="orange","D"="yellow")) {
  sapply(strings, function(x) {
    if (grepl("M1000",x) && grepl("c1",x) && !grepl("c1\\.5",x))
      return(colors["A"])
    else if (grepl("M1000",x) && (grepl("c1\\.",x) || grepl("c[2-3]",x)))
      return(colors["B"])
    else if (!grepl("M1000",x) && (grepl("c1\\.",x) || grepl("c[2-3]",x)))
      return(colors["C"])
    else if (!grepl("M1000",x) && (grepl("c1",x) && !grepl("c1\\.5",x)))
      return(colors["D"])
    else {
      print(x)
      stop ("check name input. for colSideColors short names are required.")
    }
  })
}

getShortStrings <- function(strings) {
  strings <- sapply(strings, shortSimulationName)
  strings
}

plotHeatmap.binary <- function(tableValues, tableText, 
                               zeroCol="lightblue", oneCol="green", naCol="gray" ,noteCol="red",
                               Acol="pink", Bcol="blue", Ccol="orange", Dcol="yellow",
                               ...) {
  require(gplots)
  bk <- c(0,0,1)
  colors <- colorpanel(length(bk)-1, zeroCol, oneCol)
  colnames(tableValues) <- getShortStrings(colnames(tableValues))
  rownames(tableValues) <- getShortStrings(rownames(tableValues))
  colsidecols <- getColorCode(colnames(tableValues),
                              colors=c("A"=Acol, "B"=Bcol, "C"=Ccol, "D"=Dcol))
  ylab <- substr(rownames(tableValues)[1], 1, 1)
  if (ylab=="s") ylab <- "driver selection [s_r]"
  rownames(tableValues) <- sub(sprintf("%s(.+)",ylab), "\\1", rownames(tableValues))
  heatmap.2(tableValues, breaks=bk, col=colors, key=FALSE, trace="none", dendrogram="none",
            Colv=FALSE, Rowv=FALSE, cellnote=tableText, na.col=naCol, notecol=noteCol, notecex=3.2,
            ColSideColors=colsidecols, ylab=ylab, ...)
}

# tmp1 <- getWaitingTimeTable("[ABD]", 20, getHeatmapEntry_binary, getHeatmapOverlayText_binary)
# pdf("~/tmp.pdf")
# # colCols <- sapply(colnames(tmp1[[1]]), function(x) {
# #   if
# # })
# plotHeatmap.binary(tmp1[[1]],tmp1[[2]],
#                    margins=c(8,8), cexRow=1, cexCol=1)
# legend("topleft", legend=c("Has reached cancer", "Has no reached cancer",
#                            "* Has reached mutator phenotype", "Data missing"),
#        fill=c("green", "lightblue", "red", "gray"))
# dev.off()

########################
### mutational waves ###
########################

getWavePoints <- function(obj, theclass=c(drivers=NA, superdrivers=NA,
                                          mutators=NA, passengers=NA)) {
  res <- lapply(obj, getSumOfClasses, thelist=as.list(theclass))  
  list(xxx=as.numeric(sub("g([0-9]+)", "\\1", names(res))),
       yyy=as.numeric(res))
}

getWaves <- function(obj, variants=list(drivers=0:10, superdrivers=0:5,
                                        mutators=NA, passengers=NA)) {
  res <- lapply(variants$drivers, function(i) {
    t1 <- lapply(variants$superdrivers, function(j) {
      t2 <- lapply(variants$mutators, function(l) {
        t3 <- lapply(variants$passengers, function(k) {
          print("getting wave..")
          res <- getWavePoints(obj, theclass=c(drivers=i,
                                        superdrivers=j,
                                        mutators=l,
                                        passengers=k))
          res
        })
        names(t3) <- sprintf("passengers%s", variants$passengers)
        t3
      })
      names(t2) <- sprintf("mutators%s", variants$mutators)
      t2
    })
    names(t1) <- sprintf("superdrivers%s", variants$superdrivers)
    t1
  })
  names(res) <- sprintf("drivers%s", variants$drivers)
  
  res
}

plotMutWaves <- function(
                        pdfPath="~/tmp.pdf",
                        fixedParams=c(driver=0.01, factor=2, 
                                             ratenormal=0.00000001, # 1e-8
                                             mutator=1000),
                         variants=list(drivers=0:10, superdrivers=0:5,
                                       mutators=NA, passengers=NA),
                         plotY=NULL, plotType="b", plotEveryXY=10, 
                         ylab="Number of cells", xlab="Generation",
                         plotLegend=TRUE, isCol="d2", isPch="d1", K=NULL,
                        filter2dVec_in=NA, plotGrid=TRUE, plotMidX=FALSE, 
                        midFunX=median, plotX=NULL, graphPaper=TRUE,
                        graphPaper.adaptLegend=TRUE, whichReplicate=NULL, ...) {
  
  if (!setequal(names(fixedParams), c("driver","factor","ratenormal","mutator")))
      stop("names(fixedParams) needs to setequal \"driver\", \"factor\", etc.")

  if (sum(sapply(variants, function(x) all(is.na(x)))) < 2)
    stop("at least two elements of variants needs to be NULL")
  
  if (sum(sapply(variants, function(x) all(is.na(x)))) < 1)
    stop("at least one elements of variants needs to be not NULL")
  
  if (!isCol%in%c("d1","d2"))
    stop("isCol needs to be element 'd1', 'd2'")

  if (!isPch%in%c("d1","d2"))
    stop("isPch needs to be element 'd1', 'd2'")
  
  if (isPch==isCol)
    warning("isPch and isCol should not be the same")
  
  plotY.orig <- plotY
  plotX.orig <- plotX
  
  if (!is.null(plotY)) {
#     plotY <- identity # quick fix to use log="y" parameter in plot() rather than own function
    ylab <- sprintf("%s [%s]", ylab, deparse(substitute(plotY)))
  }
  else 
    plotY <- identity
  
  ## which rdata file should be loaded
  if (is.null(whichReplicate)) {
    case <- sprintf(".*mutator%sdriver%sfactor%sratenormal%s_4mutClassesMean\\.rda", 
                    fixedParams["mutator"], fixedParams["driver"], 
                    fixedParams["factor"], fixedParams["ratenormal"])
  } else {
    case <- sprintf(".*mutator%sdriver%sfactor%sratenormal%s_4mutClasses\\.rda", 
                    fixedParams["mutator"], fixedParams["driver"], 
                    fixedParams["factor"], fixedParams["ratenormal"])
  }
  rdaname <- grep(case, list.files(intermediateDir), value=TRUE)
  rdaname <- rdaname[!grepl("final", rdaname)]

  if (length(rdaname) != 1) {
    print(case)
    stop("check case")
  }
  
  ## load rdata and get waves
  print(sprintf("loading %s", rdaname))
  obj <- loadIfNot(file.path(intermediateDir,rdaname))
  if (!is.null(whichReplicate)) {
    obj <- lapply(obj, "[[", whichReplicate)
  }
  waves <- getWaves(obj, variants)
  if (is.null(whichReplicate)) {
    savename <- sprintf("waves_%s", rdaname)
    objname <- gsub("\\.rda", "", savename)
    assign(objname, waves)
    save(list=c(objname), file=file.path(intermediateDir, "waves", savename))
  } else {
    savename <- sprintf("waves_%s_%s", whichReplicate, rdaname)
    objname <- gsub("\\.rda", "", savename)
    assign(objname, waves)
    save(list=c(objname), file=file.path(intermediateDir, "waves", savename))
  } 
  
  ## for limits of plot
  allxs <- sapply(waves, function(a) sapply(a, function(b) sapply(b, function(c) sapply(c, function(d) d$xxx))))
  allys <- sapply(waves, function(a) sapply(a, function(b) sapply(b, function(c) sapply(c, function(d) d$yyy))))
  
  xli <-range(allxs, na.rm=T)
  
  ## other Y scale if wished
#   if (plotYlog) {
      allys <- plotY(allys)
      allys[is.infinite(allys)] <- 0
#       ylab <- sprintf("%s (log scale)", ylab)
#   }
  yli <- range(allys, na.rm=T)
  yli <- c(0,yli[2]) # 0 to max actually

  pdf(pdfPath, height=16, width=21)
  
  ## plotting stuff ##
  
  if (plotLegend & !graphPaper)
    layout(rbind(1,2), heights=c(9,1)) 

  if (plotLegend & graphPaper & !graphPaper.adaptLegend) {
    parop <- par(mar=c(0, 0, 0, 0))
    layout(rbind(1,2), heights=c(7,3)) 
  } else if (plotLegend & graphPaper & graphPaper.adaptLegend) {
    parop <- par(mar=c(0, 0, 0, 0))
    layout(rbind(1,2,3), heights=c(6,0.5,1)) 
  }

  options("scipen"=1, "digits"=10)
  main <- do.call(paste, as.list(sapply(names(fixedParams), 
                                        function(x) sprintf("%s=%s", x, fixedParams[x]))))
  main <- shortSimulationName(gsub(" ", ", ", main)) # replace " " with , and make short names afterwards
  if (!is.null(plotY.orig))
    fixLog <- "y"
  else
    fixLog <- ""
  if (!is.null(plotX.orig))
    fixLog <- sprintf("%sx",fixLog)
  fixLog <- ""
  options("scipen"=100, "digits"=1)
  

if (!graphPaper) {
  plot(1,1,type="n",xlim=xli,ylim=yli,xlab=xlab,ylab=ylab,main=main, yaxt="n", ...)
  if (is.null(plotY.orig)) 
    axis(2, axTicks(2), format(axTicks(2), scientific=T))
  else
    axis(2, axTicks(2), format(axTicks(2), scientific=F)) 
} else {
  if (!graphPaper.adaptLegend)
    op <- par(mar = c(10,10,4,2) + 0.1) ## default is c(5,4,4,2) + 0.1
  else
    op <- par(mar = c(12,17,4,2) + 0.1) ## default is c(5,4,4,2) + 0.1
  if (!graphPaper.adaptLegend)
    plot(1,1,type="n",xlim=xli,ylim=yli,main="",yaxt="n",cex.axis=3, 
         ann = FALSE, tck=-0.03, padj=1.5) 
  else
    plot(1,1,type="n",xlim=xli,ylim=yli,main="",yaxt="n",cex.axis=5, 
         ann = FALSE, tck=-0.03, padj=1.5) 
  
  if (is.null(plotY.orig)) {
    if (!graphPaper.adaptLegend)
      axis(2, axTicks(2), format(axTicks(2), scientific=T),
           cex.axis=3, tck=-0.03, padj=-0.69)
    else
      axis(2, axTicks(2), format(axTicks(2), scientific=T),
           cex.axis=5, tck=-0.03, padj=-0.69)
  } else {
    if (!graphPaper.adaptLegend)
      axis(2, axTicks(2), format(axTicks(2), scientific=F), 
           cex.axis=3, tck=-0.03, padj=-0.69) 
    else
      axis(2, axTicks(2), format(axTicks(2), scientific=F), 
           cex.axis=5, tck=-0.03, padj=-0.69) 
  }
    
  if (!graphPaper.adaptLegend)
    title(ylab=ylab,line=5,xlab=xlab,cex.lab=3)
  else
    title(ylab=ylab,line=11,xlab=xlab,cex.lab=5)
  #     axis(1, at = c(0), las=2, tck=-0.045, labels=TRUE, lwd=3, cex.axis=.7)
par(op)
}
  
  if (plotGrid) 
    grid(lwd=9)
  
  ##########################
  
  ## to infer if driver, superdriver, mutator, or passengers should be analyzed in plots
  idx <- which(!sapply(variants, function(x) all(is.na(x)))) # which variants are not NA
  d1_name <- names(idx)[1]
  d2_name <- names(idx)[2]
  ## result array storing mid points, where mid is either median, or mean, etc.
  resArr <- matrix(NA, ncol=length(variants[[d2_name]]), nrow=length(variants[[d1_name]]))
  colnames(resArr) <- 0:(length(variants[[d2_name]])-1)
  rownames(resArr) <- 0:(length(variants[[d1_name]])-1)
#   dimnames(resArr) <- list(d1_name,d2_name)
  
  lg_col <- c()
  lg_pch <- c()
  lg_lty <- c()
  lg_text <- c()
  
  ## loop through all waves of mutation types
  for (i in 1:length(waves)) {
    for (j in 1:length(waves[[i]])) {
      for (k in 1:length(waves[[i]][[j]])) {
        for (l in 1:length(waves[[i]][[j]][[k]])) {
          
            # this variant whose first idx is not null determines the col
            if (idx[1]==1) {
              d1 <- i
            } else if (idx[1]==2) {
              d1 <- j
            } else if (idx[1]==3) {
              d1 <- k
            } else if (idx[1]==4) {
              d1 <- l
            } else stop("something is wrong with idx and col")
            
            # this variant whose first idx is not null determines the lty
            if (length(idx) > 1)
              d2 <- j # just take second loop
            else {
              if (idx[2]==1)
                d2 <- i
              else if (idx[2]==2)
                d2 <- j
              else if (idx[2]==3)
                d2 <- k
              else if (idx[2]==4)
                d2 <- l
              else stop("something is wrong with idx and lty")
            }
            
            d1 <- d1-1
            d2 <- d2-1
            
            if (is.null(K) || d1+d2 == K) {
              lclWaves <- waves[[i]][[j]][[k]][[l]]
              xxx <- lclWaves$xxx
              yyy <- lclWaves$yyy
              xxx <- xxx[seq(1,length(xxx), by=plotEveryXY)]
              yyy <- yyy[seq(1,length(yyy), by=plotEveryXY)]
              posEntries <- which(!yyy<=1)
              yyy <- yyy[posEntries]
              xxx <- xxx[posEntries]
              yyy <- plotY(yyy)
              
              if (isCol=="d1")
                col <- d1+1
              else
                col <- d2+1
              
              if (isPch=="d2") {
                lty <- d2+1
                pch <- d2+1
              } else {
                lty <- d1+1
                pch <- d1+1
              }
              miniVector <- c(variants[[d1_name]][d1+1], variants[[d2_name]][d2+1])
              ## next lines for when peak of wave (or left edge or so) should be computed
              tmpa <- as.character(miniVector[1])
              tmpb <- as.character(miniVector[2])
              tmpmid <- midFunX(xxx)
              resArr[tmpa, tmpb] <- tmpmid
#               if (all(miniVector==c(0,1)))
#                 browser()
              ## if filter is turned on
              if (is.na(filter2dVec_in) || any(sapply(filter2dVec_in, function(x) all(x==miniVector)))) {
#               if (is.na(filter2dVec_in) || list(miniVector)%in%filter2dVec_in) {
                lg_pch <- c(lg_pch, pch)
                lg_lty <- c(lg_lty, pch)
                lg_col <- c(lg_col, col)
                lg_text_tmp <- sprintf("(%s, %s)", miniVector[1], miniVector[2])          
                lg_text <- c(lg_text, lg_text_tmp)
                if (!graphPaper)
                  lines(xxx, yyy, col=col, lty=lty, pch=pch, type=plotType, lwd=2, ...)
                else
                  lines(xxx, yyy, col=col, lty=lty, pch=pch, type=plotType, lwd=5, cex=3, ...)
                if (plotMidX) abline(v=tmpmid, col=col, lty=lty)
              } ##
          } else print("K not null OR d1+d2 != K")
        }          
      }
    }
  }
  ## plot legend
  if (plotLegend) {
    if (!graphPaper) {
      parop <- par(mar=c(0, 0, 0, 0))
      plot.new()
      if (plotType!="b")
        lg_lty=NULL
      legend("center",
             legend = lg_text, 
             bty = "n", 
             ncol=length(unique(lg_pch)),#max(length(variants[[d1_name]]),length(variants[[d2_name]])),
             col = lg_col,
             lty = lg_lty, 
             pch = lg_pch,
             #pt.bg = lg_pt.bg,
             #pt.cex = 2
      )
      par(parop) 
    } else {
      plot.new()
      if (plotType!="b")
        lg_lty=NULL
      if (!graphPaper.adaptLegend) {
        legend("center",
               legend = lg_text, 
               bty = "n", 
               ncol=length(unique(lg_pch)),#max(length(variants[[d1_name]]),length(variants[[d2_name]])),
               col = lg_col,
               lty = lg_lty, 
               pch = lg_pch,
               #pt.bg = lg_pt.bg,
               pt.cex = 2,
               cex = 2
        )
      } else {
        print("Plot grpahPaper.adaptLegend")
        nc <- length(unique(lg_pch))
        nr <- length(unique(lg_col))
        tt <- seq(from=1,to=length(lg_col),by=nr)
        tt2 <- 1:nr
        lg_pch.cut <- lg_pch[tt]
        lg_lty.cut <- lg_lty[tt]
        lg_col.cut <- lg_col[tt2]
        ks <- gsub("\\(\\d+, (\\d+)\\)", "\\1", lg_text[1:nr])
        kr <- gsub("\\((\\d+), \\d+\\)", "\\1", lg_text[tt])
        legend("bottom", legend=paste0("ks=",ks), fill=lg_col.cut,
               pt.cex=3.5, cex=3.5, bty="n", ncol=length(ks))
        plot.new()
        legend("top", legend=paste0("kr=",kr), lty=lg_lty.cut, 
               pch=lg_pch.cut, pt.cex=3, cex=3, bty="n", ncol=length(ks))
      }
      par(parop)
    }
  }  
  dev.off() ## shut down pdf device
  resArr
}

#' @param thisIn list of miniVector (i,j), ordered according to !expected! occurence in plot
stepwiseDistances <- function(resArr, thisIn) {
  values <- unlist(sapply(thisIn, function(xy) {
    xy <- as.character(xy)
    resArr[xy[1],xy[2]]
  }))
  
  rec <- function(vals, dist=c()) {
    if (length(vals) == 1)
      return(dist)
    else {
      dist2 <- vals[2] - vals[1]
      rec(vals[-1], c(dist,dist2)) 
    }
  }
  
  rec(values)
}

#' @param filterIn list of 2d vectors (minVector)
boxplotMidX <- function(listOfStepwiseDists, ...) {  
  boxplot(listOfStepwiseDists, ...)
}

###################################
## analyze waves
###################################

getRuns <- function(yy) {
  with(rle(yy!=0), {
    ok <- values
    ends <- cumsum(lengths)[ok]
    starts <- ends - lengths[ok] + 1
    cbind(starts, ends)
  })
}

analyzeWaves <- function(waves, whichStat=c("width", "height", "peak_x-distance-superdriver", "peak_y-distance-driver", "num_driver")) {
  whichStat <- match.arg(whichStat)
  
  outerf <- function(waves, fun) {
    lapply(waves, function(driver) {
      lapply(driver, function(superdriver) {
        lapply(superdriver, function(mutator) {
          lapply(mutator, function(passenger) {
            fun(passenger)
          })
        }) 
      })
    })   
  }
  
  getYYYruns <- function(passenger) {
    xxx <- passenger$xxx
    yyy <- passenger$yyy
    mat <- getRuns(yyy)
    list(xxx=xxx, yyy=yyy, yyy_runs=mat)
  }
  
  waves.yyy_runs <- outerf(waves, getYYYruns)
  
  run.width <- function() {
    getWidth <- function(passenger) {
      xxx <- passenger$xxx
      yyy <- passenger$yyy
      mat <- passenger$yyy_runs
      if (nrow(mat) > 0) {
        ## get longest run
        lens <- apply(mat, 1, function(x) x[2]-x[1]) # end - start index
        idx <- mat[which.max(lens), ]
        ## length of max wave measured in generation (xxx)
        xxx[idx[2]] - xxx[idx[1]] 
      } else 0
    }
    res <- outerf(waves.yyy_runs, getWidth)
    unlist(res)
  }
  
  run.height <- function() {
    getHeight <- function(passenger) {
      xxx <- passenger$xxx
      yyy <- passenger$yyy
      mat <- passenger$yyy_runs
      if (nrow(mat) > 0) {
        ## get longest run
        lens <- apply(mat, 1, function(x) x[2]-x[1]) # end - start index
        idx <- mat[which.max(lens), ]
        ## length of wave measured in generation (xxx)
        max(yyy[idx[1]:idx[2]])
      } else 0
    }
    res <- outerf(waves.yyy_runs, getHeight)
    unlist(res)
  }
  
  run.num_driver <- function() {
    widths <- run.width()
    widths <- widths[widths!=0] # this is equivalent to total number of waves
    # get superdrivers superdrivers
    superdrivers <- sort(unique(gsub(".+superdrivers(\\d+).+", "\\1", names(widths))))
    superdrivers <- sprintf("superdrivers%s", superdrivers)
    sapply(superdrivers, function(sudr) sum(grepl(sudr, names(widths))))
  }
  
 
  run.peak_x <- function() {
    heights <- run.height()
    heights <- heights[heights!=0]
    # get drivers 
    drivers <- sort(as.numeric(unique(gsub("^drivers(\\d+).+", "\\1", names(heights)))))
    drivers <- sprintf("drivers%s", drivers)
    # get difference in x direction between all superdriver waves conditioned on number of drivers
    sapply(drivers, function(dr) {
      # fixed driver number
      superdriverHeights <- heights[sort(grep(sprintf("^%s\\.",dr), names(heights), val=T))]
      res <- sapply(names(superdriverHeights), function(sudr) {
        tt <- waves[[dr]][[gsub("drivers\\d+\\.(superdrivers\\d+)\\.mutators.+", "\\1", sudr)]][[1]][[1]]
        tt$xxx[which(tt$yyy==heights[[sudr]])[1]]
      })
      # res <- res[sapply(res, length) > 0]
      names(res) <- gsub("drivers\\d+\\.(superdrivers\\d+)\\.mutators.+", "\\1", names(res))
      diff(res)
    })
  }
  
  ## For a superdriver block, record difference in height of waves with increasing driver numbers 
  run.peak_y <- function() {
    heights <- run.height()
    heights <- heights[heights!=0]
    # get superdrivers 
    superdrivers <- sort(unique(gsub(".+superdrivers(\\d+).+", "\\1", names(heights))))
    superdrivers <- sprintf("superdrivers%s", superdrivers)
    # get difference in peaks of drivers for every superdriver block
    res <- sapply(superdrivers, function(sudr) diff(heights[sort(grep(sudr, names(heights), val=T))]))
    res <- sapply(res, function(x) {
      names(x) <- gsub("^(drivers\\d+)\\.super.+", "\\1", names(x))
      x
    })
    res
  }
  
  switch(whichStat,
         "width"=run.width(),
         "height"=run.height(),
         "num_driver"=run.num_driver(),
         "peak_x-distance-superdriver"=run.peak_x(),
         "peak_y-distance-driver"=run.peak_y())
}


plotWaveStats <- function(cases="^waves_B_mutator1000driver0.01factor(1\\.1|1\\.3|1\\.5|2\\.6|2\\.8|3)ratenormal0.00000001_4mutClassesMean.rda") {
  files <- list.files(file.path(intermediateDir, "waves"), pattern=cases, full.names=T)
  for (file in files) load(file)
  params <- gsub("\\.rda", "", basename(files))
  stats <- c("width", "height", "peak_x-distance-superdriver", "peak_y-distance-driver", "num_driver")
  res <- lapply(stats, function(stat) {
    tt <- lapply(params, function(waves) {
      xx <- get(waves)
      analyzeWaves(xx, whichStat=stat)
    })
    names(tt) <- params 
    tt
  })
  names(res) <- stats
  res

  plotNumDrivers <- function(res) {
    ## Barplot where x axis is c and y axis is number of driver waves
    ## Every x-entry has multiple bars (one per superdriver)
    dfs <- lapply(names(res), function(nam) {
      fac <- gsub(".+factor(.+)rate.+", "\\1", nam)
      data.frame("NumDrivers"=res[[nam]], "NumSuperdrivers"=gsub("superdrivers", "", names(res[[nam]])), "Factor"=fac)
    })
    df <- Reduce(rbind, dfs)    
    library(ggplot2)
    gp <- ggplot(df, aes(x=Factor, y=NumDrivers, fill=NumSuperdrivers)) 
    gp <- gp + geom_bar(stat="identity",position="dodge")
    gp <- gp + xlab("Superdriver selection")  + ylab("Number of drivers") + scale_fill_discrete(name="# superdrivers")
    gp <- gp + theme_minimal(base_size=21)
    #gp <- gp +  scale_fill_discrete(name="Gender", breaks=c(1, 2), labels=c("Male", "Female")) 
    ggsave(file.path(plotDir, "waveanalytics_numdrivers.pdf"), gp) 
  }

  plotWidths <- function(res, min.generation=1100) {
    dfs <- lapply(names(res), function(nam) {
      fac <- gsub(".+factor(.+)rate.+", "\\1", nam)
      vec <- res[[nam]]
      vec <- vec[vec!=0]
      data.frame("Width"=vec, "Factor"=fac)
    })
    df <- Reduce(rbind, dfs)
    df <- df[df$Width > min.generation, , drop=F]
    require(ggplot2)
    gp <- ggplot(df, aes(x=Factor, y=Width)) + geom_boxplot() + geom_jitter(width=0.1)
    gp <- gp + xlab("Superdriver selection") + ylab("Width of waves") 
    gp <- gp + theme_minimal(base_size=21)
    ggsave(file.path(plotDir, "waveanalytics_widths.pdf"), gp)

  }
   
  plotYDistances <- function(res) {
    dfs <- lapply(names(res), function(nam) {
      fac <- gsub(".+factor(.+)rate.+", "\\1", nam)
      dfs_sub <- lapply(names(res[[nam]]), function(sudr) {
        sudr_res <- res[[nam]][[sudr]]
        data.frame(NumSuperdriver=gsub("superdrivers", "\\1", sudr), PeakDifference=sudr_res, Factor=fac)
      })
      Reduce(rbind, dfs_sub)
    })
    df <- Reduce(rbind, dfs)
    require(ggplot2)
    gp <- ggplot(df, aes(x=Factor, y=PeakDifference, fill=NumSuperdriver)) + geom_boxplot(outlier.shape=NA) 
    #gp <- gp + scale_y_log10() 
    gp <- gp + coord_cartesian(ylim = c(-25000000, 50000000))
    gp <- gp + xlab("Superdriver selection") + ylab("Peak difference") 
    gp <- gp + theme_bw(base_size=21) + scale_fill_discrete(name="# superdrivers")
    gp <- gp + scale_y_continuous(labels = function(x) format(x, scientific = TRUE))
    gp <- gp + theme(axis.text.y = element_text(angle = 90, hjust = 1))
    ggsave(file.path(plotDir, "waveanalytics_difference_y.pdf"), gp)

  }

  plotXDistances <- function(res, max.drivers=10) {
    dfs <- lapply(names(res), function(nam) {
      fac <- gsub(".+factor(.+)rate.+", "\\1", nam)
      res_oneRun <- res[[nam]]
      dfs_sub <- lapply(names(res_oneRun), function(dr) {
        res_dr <- res_oneRun[[dr]]
        if (length(res_dr) > 0)
          data.frame(NumDriver=gsub("drivers", "", dr), GenDifference=res_dr, Factor=fac)
        else data.frame(NumDriver=gsub("drivers", "", dr), GenDifference="EMPTY", Factor=fac)
      })
      Reduce(rbind, dfs_sub)
    })
    df <- Reduce(rbind, dfs)
    df <- subset(df, GenDifference != "EMPTY" & GenDifference != 0)
    df$GenDifference <- as.numeric(df$GenDifference)
    df$NumDriver <- as.numeric(df$NumDriver)
    df <- df[df$NumDriver <= max.drivers, , drop=F]
    df$NumDriver <- factor(df$NumDriver)
    require(ggplot2)
    gp <- ggplot(df, aes(x=Factor, y=GenDifference, fill=NumDriver)) + geom_boxplot(outlier.shape=NA)
    gp <- gp + xlab("Superdriver selection") + ylab("Peak distance")
    gp <- gp + theme_bw(base_size=21) + scale_fill_discrete(name="# drivers")
    ggsave(file.path(plotDir, "waveanalytics_difference_x.pdf"), gp)
  }

  plotXDistances(res[["peak_x-distance-superdriver"]])
  plotYDistances(res[["peak_y-distance-driver"]])
  plotWidths(res$width)
  plotNumDrivers(res$num_driver)
  
}

### More sophisticated stats
getQuadraticParams <- function(x, y) {

  res.lm <- lm(y ~ x + I(x^2), data=data.frame(y=y, x=x))
  res.lm.summary <- summary(res.lm)
  c <- res.lm.summary$coefficients["I(x^2)", "Estimate"]
  b <- res.lm.summary$coefficients["x", "Estimate"]
  a <- res.lm.summary$coefficients["(Intercept)", "Estimate"]

  curvature <- sign(c) * sqrt(abs(c))
  location_temp <- (b/2/curvature)
  location <- -1 * sign(c) * location_temp/curvature
  height <- a - sign(c) * location_temp^2
  
  list(params=c(curvature=curvature, location=location, height=height),
       coefficients=c(inter=a, x=b, x2=c),
       model=list(res=res.lm, summary=res.lm.summary),
       input=list(xxx=x, yyy=y))

}

getQuadraticDistribution <- function(case, replicates=1:50, drivers=0:10, superdrivers=0:5, noReplicate=F) {
  
  if (!noReplicate)
    cases <- paste0("waves_", replicates, "_", case, ".rda")
  else {
    replicates <- 1
    cases <- paste0("waves_", case, "Mean.rda")
  }
  
  waves <- lapply(cases, function(x) {
    x <- file.path(intermediateDir, "waves", x)
    if (file.exists(x)) {
      load(x, envir=.GlobalEnv)
      get(gsub("\\.rda", "", basename(x)))
    } else NULL
  })
  waves <- waves[!is.null(waves)]
  
  res <- lapply(drivers, function(dr) {
    res.superdrivers <- lapply(superdrivers, function(su) {
      res.waves <- lapply(replicates, function(rep) {
        wave <- waves[[rep]][[paste0("drivers",dr)]][[paste0("superdrivers",su)]][["mutatorsNA"]][["passengersNA"]]
        if (all(wave$yyy==0)) {
          return(NULL)
        } else {
          mat <- getRuns(wave$yyy)
          lens <- apply(mat, 1, function(x) x[2]-x[1]) # end - start index
          idx <- mat[which.max(lens)[1],]
          num_gen <- (idx[2]-idx[1]) * 10
          if (! num_gen > 100) {
            return(NULL)
          } else {
            xxx <- wave$xxx[idx[1]:idx[2]]
            yyy <- wave$yyy[idx[1]:idx[2]]
            tryCatch(getQuadraticParams(xxx, yyy), error=function(e) browser())
          }
        }
      })
      names(res.waves) <- paste0("replicate", replicates)
      res.waves
    })
    names(res.superdrivers) <- paste0("superdrivers", superdrivers)
    res.superdrivers
  })
  names(res) <- paste0("drivers", drivers)
  
  res.compact <- lapply(res, function(dr) lapply(dr, function(su) Reduce(rbind, lapply(su, "[[", "params"))))
  
  res.compact
  
}

plotQDistributions <- function(cases, plot.drivers=0:2, plot.superdrivers=0:3, plotAsterisk=T, ...) {
  res.all <- lapply(cases, getQuadraticDistribution, ...)
  names(res.all) <- cases
  
  res.all_mean <- lapply(cases, getQuadraticDistribution, noReplicate=T, ...)
  names(res.all_mean) <-cases
  
  ## We make three multi plot figure, where x axis are drivers and y are superdrivers
  ## Within everyy plot, x are superdriver selection and y are either curvature, location, or height
  
  createTable <- function(res=res.all, col=c("curvature", "location", "height")) {
    col <- match.arg(col)
    
    tt <- lapply(names(res), function(fac.name) {
      fac <- res[[fac.name]]
      res.fac <- lapply(names(fac[paste0("drivers",plot.drivers)]), function(dr.name) {
        dr <- fac[[dr.name]]
        res.dr <- lapply(names(dr[paste0("superdrivers",plot.superdrivers)]), function(su.name) {
          su <- dr[[su.name]]
          if (!is.null(dim(su)))
            distr <- su[, col]
          else distr <- su[col]
          data.frame("Factor"=gsub("^.+(factor[123568\\.]+)ratenormal.+$", "\\1", fac.name),
                     "Driver"=as.character(dr.name),
                     "Superdriver"=as.character(su.name),
                     "Value"=distr,
                     "Type"=col)
        })
        names(res.dr) <- names(dr[plot.superdrivers])
        Reduce(rbind, res.dr)
      })
      names(res.fac) <- names(fac[plot.drivers])
      Reduce(rbind, res.fac)
    })
    names(tt) <- names(res)
    Reduce(rbind, tt)
  }
  
  tab.curvature <- createTable(res.all, "curvature")
  tab.location <- createTable(res.all, "location")
  tab.height <- createTable(res.all, "height")
  
  tab.curvature_mean <- createTable(res.all_mean, "curvature")
  tab.location_mean <- createTable(res.all_mean, "location")
  tab.height_mean <- createTable(res.all_mean, "height")
  
  doPlot <- function(tab, name, tab_mean) {
    require(ggplot2)
    if (name=="Location") {
      tab <- subset(tab, Driver!="drivers0")
      tab_mean <- subset(tab_mean, Driver!="drivers0")
    }
    ## Change formatting a little
    tab$Factor <- as.factor(as.character(gsub("factor", "", tab$Factor)))
    tab$Driver <- gsub("drivers", "D", tab$Driver)
    tab$Superdriver <- gsub("superdrivers", "S", tab$Superdriver)
    tab_mean$Factor <- as.factor(as.character(gsub("factor", "", tab_mean$Factor)))
    tab_mean$Driver <- gsub("drivers", "D", tab_mean$Driver)
    tab_mean$Superdriver <- gsub("superdrivers", "S", tab_mean$Superdriver)
    ## Start plotting
    gp <- ggplot(tab, aes(x=Factor, y=Value)) + geom_boxplot(outlier.shape=NA)
    gp <- gp + xlab("Superdriver selection") + ylab(name)
    gp <- gp + theme_bw(base_size=21)
    gp <- gp + facet_grid(Driver ~ Superdriver, scales="free")
    gp <- gp + theme(axis.text.x=element_text(angle=90, hjust=1))
    if (name=="Height")
      gp <- gp + scale_y_log10(labels = function(x) format(x, scientific = TRUE))
    # if (name=="Height" | name=="Location")
    #   gp <- gp + scale_y_continuous(labels = function(x) format(x, scientific = TRUE))
    if (plotAsterisk)
      gp <- gp + geom_text(label="*", data=tab_mean, size=7, col="blue")
    gp <- gp + theme(plot.title = element_text(hjust = 0.5)) #+ ggtitle(name) 
    ggsave(file.path(plotDir, sprintf("wavestats_%s.pdf", name)), gp, width=12, height=8, dpi=300)
  }
  
  doPlot(tab.curvature, "Curvature", tab_mean=tab.curvature_mean)
  doPlot(tab.location, "Location", tab_mean=tab.location_mean)
  doPlot(tab.height, "Height", tab_mean=tab.height_mean)
  
  
}

###################################
## attribute vs attribute over time
###################################

## unfinished
attrVSattrInTime <- function(fixedParams=c(driver=0.01, mutator=1000, factor=2, ratenormal=0.00000001),
                             generations=c(1, seq(1000,4500,by=1000), 4500)) {
  case <- sprintf(".*mutator%sdriver%sfactor%sratenormal%s_oneClonesRes\\.rda", 
                  fixedParams["mutator"], fixedParams["driver"], 
                  fixedParams["factor"], fixedParams["ratenormal"])
  rdaname <- grep(case, list.files(intermediateDir), value=TRUE)
  rdaname <- rdaname[!grepl("final", rdaname)]
  
  print(sprintf("loading %s", rdaname))
  obj <- loadIfNot(file.path(intermediateDir,rdaname))
  
  obj$wav 
}

###################
## below depracated
###################

plotList_points <- function(results, xParam, onPlot, ...)
  plotList(results, xParam, onPlot, mean, points, ...)

plotFor_2onPlot <- function(results)

#' generic functions
plotForFixedParams <- function(results, funFixedParams, funOnPlotParams, plotFun, centerFun, ...) {
  
  toPlotResults <-  2
  plotFun(toPlotResults, ...)
}

plotForFixedParams_linesANDsd <- function(results, fixedParams, onPlotParams, centerFun, ...) {
  
  
  plotForFixedParams(results, fixedParams, function()1)
}

#' generic function.
#' for all simulation results that match a pattern (e.g. for all folder of case A), vary one parameter
#' and fix rest (two)
#' function loads corresponding rdata files if not loaded yet and applys f on them
#' @param pattern string pattern specifying which simulation results to use among all
#' @param f function to apply (usually producing one combined plot) to the set of simulation results,
#' where only one parameter varys (usually x-axis)
plotFun_twoFixed_forAllCase <- function(pattern, f, ...) {
  ## get all *t_* file names that specify pattern normally in intermediate folder
  files <- grep(pattern, list.files(intermediateDir), value=TRUE) 
  #   files <- sub("(.+)\\.rda", "\\1", files) # get real names
  ## load files into workspace
  objects <- lapply(sapply(files, function(x) file.path(intermediateDir,x)), loadIfNot)
  names(objects) <- sub("(.+)\\.rda", "\\1", files) 
  
  combTable <- getCombTable(names(objects))
  
  performForCombs(getCombTable(names(objects)), f, ..., objects=objects)
}

#' plot waiting times, one plot for each parameter fixing for all possible combinations
#' of fixing one parameter and varying the rest. each possible parameter is varied.
plotWaitingTimes_twoFixed_forAllCase <- function(pattern)
  plotFun_twoFixed_forAllCase(pattern, plotWaitingTimes_twoFixed)

plotNclasses_twoFixed_forAllCase <- function(pattern)
  plotFun_twoFixed_forAllCase(pattern, plotNclasses_twoFixed)

# plotWaitingTimes_twoFixed_forAllCase <- function(pattern) {
#   # get all *t_* file names that specify pattern normally in intermediate folder
#   files <- grep(pattern, list.files(intermediateDir), value=TRUE) 
# #   files <- sub("(.+)\\.rda", "\\1", files) # get real names
#   # load files into workspace
#   objects <- lapply(sapply(files, function(x) file.path(intermediateDir,x)), loadIfNot)
#   names(objects) <- sub("(.+)\\.rda", "\\1", files) 
#   
#   combTable <- getCombTable(names(objects))
#   
#   perform
#   
#   for (i in c("mutator", "driver", "factor")) { #for varying mutator, driver, and factor 
#     lclCombTable <- combTable
#     lclCombTable[,i] <- ".*"
#     lclCombTable <- unique(lclCombTable)
#     apply(lclCombTable, 1, function(x) {
# #       tmpObj <- objects[findNames(names(objects), # get objects that has that mutators, drivers, factors
# #                                   mutator=x["mutator"],
# #                                   driver=x["driver"],
# #                                   factor=x["factor"])]
#       plotWaitingTimes_twoFixed(objects, 
#                                 mutator=x["mutator"],
#                                 driver=x["driver"],
#                                 factor=x["factor"])
#     })
#   }
# }

#' generic function to apply a plotting/analysis function for varying one 
#' specific parameter and fixing the rest. usually produes one single plot
#' @param results list of pre-filtered simulation results. list corresponds
#' to all possible combinations of fixed two values and specified variable
#' parameter
#' @param mutator numeric or ".*" if should be varied
#' @param driver numeric or ".*" if should be varied
#' @param factor numeric or ".*" if should be varied
#' NOTE: only one of mutator, driver, factor, can be ".*"
#' @param f plotting functions applied on results
#' @param ... additional parameters to f
plot_twoFixed <- function(results, mutator, driver, factor, ratenormal, f, ...) {
  #whichVaries1, value1, whichVaries2, value2) {
  
  if (! ".*"%in%c(mutator, driver, factor, ratenormal))
    stop("'.*' needs to be specified in one of the parameters")
  
  if (mutator==".*")
    whichVaries <- "mutator"
  else if (driver==".*")
    whichVaries <- "driver"
  else if (ratenormal==".*")
    whichVaries <- "factor"
  else
    whichVaries <- "ratenormal"
  
  #   if (! whichVaries%in%c("mutator", "driver", "factor") && 
  #         whichVaries2%in%c("mutator", "driver", "factor")))
  #     stop("'whichVaries1' or 2 not %in% 'driver', 'mutator', 'factor'")
  
  ## filter out null entries where no waiting time could be obtained
  results <- Filter(function(x) !is.null(x), results) # filter out null results
  
  # usually one parameter is ".*"
  thenames <- findNames(names(results),
                        mutator=mutator, driver=driver, factor=factor, ratenormal=ratenormal) 
  results <- results[thenames] # get only results with specified values
  
  theparams <- getParamsFromStrings(thenames) # get parameter values of matched simulation results
  # get values of varying parameter as numeric values and use as names for matched simulations results
  # in order to be plottable in boxplots etc.
  theVaryParams <- sapply(theparams, function(x) as.numeric(x[whichVaries])) 
  names(results) <- theVaryParams
  
  # order from low to high values
  theorder <- order(theVaryParams, decreasing=FALSE)
  f(results[theorder], xlab=whichVaries,
    main=sprintf("M=%s, s=%s, c=%s, u=%s", mutator, driver, factor, ratenormal), ...)
}

#' for a set of simulation results where exactly one parameter is varied, plot
#' number of mutational classes
plotNClasses_twoFixed <- function(mutClassRes, ..., mutator, driver, factor, ratenormal) {
  plot_twoFixed(mutClassRes, mutator=mutator, driver=driver, factor=factor, ratenormal=ratenormal,
                f=plotNclasses_box, ...)
}

#' depracated don't use!
#' plot number of classes for all simulation results
plotNClasses_twoFixed_old <- function(results, mutator, driver, factor, ...) {
  if (! ".*"%in%c(mutator, driver, factor))
    stop("'.*' needs to be specified in one of the parameters")
  
  if (length(which(c(mutator, driver, factor)%in%".*")) > 1)
    stop("'.*' can only be specified once")
  
  if (mutator==".*")
    whichVaries <- "mutator"
  else if (driver==".*")
    whichVaries <- "driver"
  else 
    whichVaries <- "factor"
  
  #   if (! whichVaries%in%c("mutator", "driver", "factor") && 
  #         whichVaries2%in%c("mutator", "driver", "factor")))
  #     stop("'whichVaries1' or 2 not %in% 'driver', 'mutator', 'factor'")
  
  ## filter out null entries where no waiting time could be obtained
  results <- Filter(function(x) !is.null(x), results) # filter out null results
  
  thenames <- findNames(names(results), mutator=mutator, driver=driver, factor=factor)
  results <- results[thenames] # get only results with specified values
  
  theparams <- getParamsFromStrings(thenames)
  theVaryParams <- sapply(theparams, function(x) as.numeric(x[whichVaries]))  
  names(results) <- theVaryParams
  
  theorder <- order(theVaryParams, decreasing=FALSE)
  
  res2plot <- results[theorder]
  
  plotNclasses_box(res2plot, 
                   xlab=whichVaries,
                   main=sprintf("M=%s, s=%s, c=%s", mutator, driver, factor))
}

#' for a set of simulation results where exactly one parameter is varied, plot
#' waiting times
plotWaitingTimes_twoFixed <- function(waitingTimesRes, ..., mutator, driver, factor, ratenormal) {
  plot_twoFixed(waitingTimesRes, mutator=mutator, driver=driver, factor=factor, ratenormal=ratenormal,
                f=plotWaitingTimes_box, ...)
}

#' depracated don't use!
plotWaitingTimes_twoFixed_old <- function(waitingTimesRes, ..., mutator, driver, factor) {
  #whichVaries1, value1, whichVaries2, value2) {
  
  if (! ".*"%in%c(mutator, driver, factor))
    stop("'.*' needs to be specified in one of the parameters")
  
  if (mutator==".*")
    whichVaries <- "mutator"
  else if (driver==".*")
    whichVaries <- "driver"
  else 
    whichVaries <- "factor"
  
  #   if (! whichVaries%in%c("mutator", "driver", "factor") && 
  #         whichVaries2%in%c("mutator", "driver", "factor")))
  #     stop("'whichVaries1' or 2 not %in% 'driver', 'mutator', 'factor'")
  
  ## filter out null entries where no waiting time could be obtained
  waitingTimesRes <- Filter(function(x) !is.null(x), waitingTimesRes) # filter out null results
  
  thenames <- findNames(names(waitingTimesRes),
                        mutator=mutator, driver=driver, factor=factor) # usually one parameter is ".*"
  waitingTimesRes <- waitingTimesRes[thenames] # get only results with specified values
  
  theparams <- getParamsFromStrings(thenames) # get parameter values of matched simulation results
  # get values of varying parameter as numeric values and use as names for matched simulations results
  # in order to be plottable in boxplots etc.
  theVaryParams <- sapply(theparams, function(x) as.numeric(x[whichVaries])) 
  names(waitingTimesRes) <- theVaryParams
  
  # order from low to high values
  theorder <- order(theVaryParams, decreasing=FALSE)
  plotWaitingTimes_box(waitingTimesRes[theorder], 
                       xlab=whichVaries,
                       main=sprintf("M=%s, s=%s, c=%s", mutator, driver, factor))
}

plotWaitingTimes_box <- function(results, ...) {
  if (length(results) == 0)
    warning("results is of length 0")
  
  thelist <- lapply(results, "[[", "repGens")
  
  tryCatch(
    boxplot(thelist, ...),
    error = function(e) e,
    finally = NULL
  )
  
}


plotNclasses_box <- function(results, ...) {
  if (length(results) == 0)
    warning("results is of length 0")
  
  thelist <- lapply(results, function(x) length(which(x > 0)))
  
  tryCatch(
    boxplot(thelist, ...),
    error = function(e) e,
    finally = NULL
  )
  
}


###################################################################
### Diversity, Friday Dec 16, 2016
###
### Patrick Grossmann, patrick@jimmy.harvard.edu
###
####################################################################

diversityIndex <- function(counts, method=c("simpson")) {
  method <- match.arg(method)
  
  simpson <- function() {
    n <- sum(counts)
    1-sum((counts*(counts-1))/(n*(n-1)))
  }
  
  switch(method,
         simpson=simpson())
}

diversityIndex.stats <- function(countss, method=c("simpson"), parallel=goParallel, ncores=nCores) {
  if (!goParallel)
    res <- lapply(countss, diversityIndex, method=method)
  else {
    print("goParallel diversity!")
    require(parallel)
    res <- mclapply(countss, diversityIndex, method=method, mc.cores=ncores)
  }
  stats <- unlist(res)
  list(mean=mean(stats), sd=sd(stats), stats=stats)
}

getFinalPops <- function(cases, ncores=nCores) {
  files <- list.files(rdataDir, pattern=sprintf("finalPops_%s\\.rda", cases), full.names=T)
  
  print("Load cases...")
  for (file in files) {
    print(file)
    obj <- gsub("\\.rda", "", file)
    if (!obj%in%ls(.GlobalEnv))
      load(file, envir=.GlobalEnv)
    else print("(Skip loading)")
  }
  print("Done loading.")
  finalPops <- lapply(basename(files), function(x) get(gsub("\\.rda", "", x), envir=.GlobalEnv))
  names(finalPops) <- gsub("\\.rda", "", basename(files))
  finalPops
}

getDiversity <- function(cases, xaxis="factor", ...) {
  print("\n#--- Get cases ---#")
  finalPops <- getFinalPops(cases=cases, ...)
  finalPops_replicateCounts <- lapply(finalPops, function(case) lapply(case, "[[", "count"))
  print("\n#--- Get diversity ---#")
  finalPops_diversityStats <- lapply(finalPops_replicateCounts, diversityIndex.stats, ...)
  means <- sapply(finalPops_diversityStats, "[[", "mean")
  sd <- sapply(finalPops_diversityStats, "[[", "sd")
  xx <- as.numeric(gsub(sprintf("^.+%s([[:digit:]\\.]+).+$", xaxis), "\\1",names(finalPops_diversityStats)))
  idx <- order(xx)
  list(xx=xx[idx], means=means[idx], sds=sd[idx])
}

plotDiversity <- function() {

  res_0.005 <- getDiversity("B_mutator1000driver0.005factor.*ratenormal0.00000001")
  res_0.01 <- getDiversity("B_mutator1000driver0.01factor.*ratenormal0.00000001")
  res_0.02 <- getDiversity("B_mutator1000driver0.02factor.*ratenormal0.00000001")
  res_0.025 <- getDiversity("B_mutator1000driver0.025factor.*ratenormal0.00000001")
  res_0.03 <- getDiversity("B_mutator1000driver0.03factor.*ratenormal0.00000001")
  res_0.04 <- getDiversity("B_mutator1000driver0.04factor.*ratenormal0.00000001")
  res_0.05 <- getDiversity("B_mutator1000driver0.05factor.*ratenormal0.00000001")

  pdf("analysis/scripts/test.pdf")

  plot(res_0.005$xx, res_0.005$means, type="b", pch=20, ylim=c(0,1))
  points(res_0.01$xx, res_0.01$means, col="green", pch=20, type="b")
  points(res_0.02$xx, res_0.02$means, col="red", pch=20, type="b")
  points(res_0.025$xx, res_0.025$means, col="brown", pch=20, type="b")
  points(res_0.03$xx, res_0.03$means, col="orange", pch=20, type="b")
  points(res_0.04$xx, res_0.04$means, col="blue", pch=20, type="b")
  points(res_0.05$xx, res_0.05$means, col="yellow", pch=20, type="b")

  legend("bottomleft", legend=paste0("s = ", c(0.005, 0.01, 0.02, 0.025, 0.03, 0.04, 0.05)), 
         fill=c("black", "green", "red", "brown", "orange", "blue", "yellow"))

  dev.off()
}
