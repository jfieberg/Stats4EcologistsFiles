#' # Maximum Likelihood Example:  Slug 'count' data set

#' Load libraries
#+warning=FALSE, message=FALSE
  library(Data4Ecologists)

#' Read in data (and write them to output data file)
  data(slugs) 
  
#' ## Negative Binomial distribution
  
#' Log likelihood for Negative Binomial Distribution
  minus.logL.s<-function(pars, slugs){
    mu<-pars[1]
    theta<-pars[2]
    -sum(dnbinom(slugs, mu=mu, size=theta, log=TRUE))  
  }
  
#' ## Finding the MLE   
#' 
#' 
#' Optim can be used to find the parameters that minimize any function.  It takes
#' 
#' - a vector of "starting" parameters (first guesses) as its first argument. 
#' - a function to minimize as its second argument (the first argument of this function must correspond to the vector
#'   of parameters that we want to estimate)
#' - you can send other "stuff" (like the data, here)
#' - You can choose from many different methods (here, we will use BFGS)
#' - We will specify that we want it to also return the Hessian (see lecture notes and code that follows)
#'
#' Note: we may get some warnings here about NaN's being produced.  The parameter lambda of the poisson 
#' distribution has to be positive.  If the method "tries" out negative values when searching for the optim,
#' it will result in NaN's and thus, the warnings.
  mle.fit<-optim(par=c(2, 10), fn=minus.logL.s,  method="BFGS", hessian=T, slugs=slugs$slugs)  
  mle.fit 
  
#' Notes on the values returned:
#' 
#' - convergence = 0 suggests we successfully found a minimum
#' - value = the -logL value at the minimum (so -value is the logL at the maximum)
  
#' We can get the variance by calculating the inverse of the Hessian using "solve"  
  (varcovmat<-solve(mle.fit$hessian)) #var cov matrix of mu^ and theta^
 
#' SE's
  sqrt(diag(varcovmat))  
  
#' ## Footer
sessionInfo()  