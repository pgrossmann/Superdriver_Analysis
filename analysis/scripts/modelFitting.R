#################################################################
###
### this file contains functions and variables to test models 
###
#################################################################

#' perform nls regression, and nls2 with brute-force if nls fails
#' @param stringFormula string of formula
#' @param dat data.frame[x,y]
#' @param sta list of start values
#' @return fit
tryNlsFailNls2 <- function(stringFormula, dat, sta) {
  require(nls2)
  formu <- as.formula(stringFormula)
  tryCatch(
    nls(formu, data=dat, start=sta),
    error = function(e)
      nls2(formu, data=dat, start=sta, algorithm="brute-force")
  )  
}

#' get some models to test
#' @param x vector of x-values
#' @param y vector of y-values
models <- function(x,y) {
  dat <- data.frame(x,y) #dat data.frame[x,y] (not matrix)
  list(lm(y~x, data = dat), 
       lm(y~I(1/x), data=dat),
       lm(y ~ log(x), data = dat),
#        lm(y ~ exp(x), data=dat),
       lm(y ~ I(-exp(x)), data=dat),
       tryNlsFailNls2("y ~ I(1/x*a) * b", dat=dat, sta=list(a=1,b=1)),
       tryNlsFailNls2("y ~ I(1/x*a) + b*x", dat=dat, sta=list(a=1,b=1)),
#        nls(y ~ I(1/x*a) + b*x, data = dat, start = list(a = 1, b = 1), model=TRUE), 
       tryNlsFailNls2("y ~ I(a + b*log(x))", dat=dat, sta=setNames(coef(lm(y ~ log(x), data=dat)), c("a", "b"))),
#        nls(y ~ I(a + b*log(x)), data=dat, start = setNames(coef(lm(y ~ log(x), data=dat)), c("a", "b"))),
       tryNlsFailNls2("y ~I (exp(1)^(a + b * x))", dat=dat, sta = list(a=1,b=1)),
#        nls(y ~I (exp(1)^(a + b * x)), data=dat, start = list(a=1,b=1)),
       tryNlsFailNls2("y ~ I(a*exp(x)/b*x)", dat=dat, sta=list(a=1,b=1)),
#        nls2(y ~ I(a*exp(x)/b*x), data=dat, start=list(a=1,b=1), algorithm="brute-force"),
       tryNlsFailNls2("y ~ I(a*exp(x)^2/x*b)", dat=dat, sta=list(a=15,b=1e-5)), # tried out manually
       tryNlsFailNls2("y ~ I(a*exp(x)/x*b)", dat=dat, sta=list(a=15,b=1e-5)), # tried out manually
       tryNlsFailNls2("y ~ I(a*exp(x)^c/x*b)", dat=dat, sta=list(a=15,b=1e-5,c=3)), # tried out manually
#        nls2(y ~ I(a*exp(x)^2/x*a), data=dat, start=list(a=15,b=1e-5), algorithm="brute-force"), # tried out manually
       # next, nikos model
       tryNlsFailNls2("y ~ I( a*( log(x/b)^2 / (x*c) ))", dat=dat, sta=list(a=15,b=1e-7*110,c=0.01*log(1e6*1e9))),
#        nls2(y ~ I( a*( log(x/b)^2 / (x*c) )), data=dat,
#            start=list(a=15,b=1e-7*110,c=0.01*log(1e6*1e9)), algorithm="brute-force"), # nikos model
       tryNlsFailNls2("y ~ I(1/x*a)+b", dat=dat, sta = list(a=1,b=1)))
#        nls(y ~ I(1/x*a)+b, data=dat, start = list(a=1,b=1)))       
}

getFits <- models
          
getFits_symbolicRegression <- function(x,y) {
  require(rgp)
  dat <- data.frame(x=x, y=y)
  symbolicRegression(y ~ x,                            
                     data=dat, functionSet=mathFunctionSet,
                     stopCondition=makeStepsStopCondition(1000))
}

getBestFitFunction_symbolicRegression <- function(fits) {
  fits$population[[which.min(sapply(fits$population, fits$fitnessFunction))]]
}

predict_symbolicRegression <- function(fits) {
  
}

nikoModel <- function(t_k,k,s,u,d=110,Ninit=1e6,Nfin=1e9) {
  x <- s
  y <- t_k
  dat <- data.frame(x,y)
  nls(y ~ k*( log(s/(u*d))^2 / (s*log(Ninit*Nfin)) ) )
}

#' The Title
#' Title2
#' @param fits list of model fits
#' @return table (length(fits) x 3), with columns AICc, AIC, and the model
aicTable <- function(fits) {
  require(AICcmodavg)
  require(plyr)
  require(stringr)
  ldply(fits, function(mod){ data.frame(AICc = AICc(mod), AIC = AIC(mod), model = deparse(formula(mod))) })
}

#' @param fits list of model fits
#' @return returns best fit according to AIC (min.)
getBestFit <- function(fits) {
  tab <- aicTable(fits)
  i <- which.min(tab[,"AIC"])
  fits[[i]]
}

fitAndGetBest <- function(x,y) getBestFit(getFits(x,y))

# x <- c(0.001,0.005,0.01,0.05,0.1)
# y <- c(4000,3000,2000,1000,500)
# tmp1 <- models(x,y)
# library(AICcmodavg); library(plyr); library(stringr)
# resultTable <- ldply(tmp1, function(mod){ data.frame(AICc = AICc(mod), AIC = AIC(mod), model = deparse(formula(mod))) })
# print(resultTable)
# print(getBestFit(tmp1))
# print(predict(tmp1[[1]], data.frame(x=1:100)))
# stop("bas")

#' @param s vector of x-values (driver advantages)
nikoLine <- function(k,s,u,d=110,Ninit=1e6,Nfin=1e9) {
  k*( log(s/(u*d))^2 / (s*log(Ninit*Nfin)) )
}


sse <- function(strinFormulas, x, y) #sort(sapply(rhs, function(rhs, x, y, verbose = TRUE) {
  lapply(rhs, function(rhs, x, y, verbose = TRUE) {
    fo <- as.formula(paste("y", rhs, sep = "~"))
    nms <- setdiff(all.vars(fo), c("x", "y"))
    start <- as.list(setNames(rep(1, length(nms)), nms))
    fm <- nls(fo, data.frame(x, y), start = start)
    if (verbose) { print(fm); cat("---\n") }
    #   deviance(fm)
    fm
}, x = x, y = y)

# rhs <- c("a*x+b", "a*x*x+b*x+c", "log(x)") #"I(1/x*a)+b", "I( a*( log(x/b)^2 / (x*c) ))"

# tmp2 <- sse(rhs, x, y)
# ldply(tmp2, function(mod){ data.frame(AICc = AICc(mod), AIC = AIC(mod), model = deparse(formula(mod))) })
