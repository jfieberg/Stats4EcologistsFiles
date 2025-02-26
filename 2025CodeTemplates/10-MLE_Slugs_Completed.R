#' # Maximum Likelihood Example:  Slug 'count' data set

#' Remove anything in memory
  rm(list = ls())

#' Load libraries
#+warning=FALSE, message=FALSE
  library(Data4Ecologists)
  
#' Read in data 
  data(slugs)
   
#' Lets start off by assuming the distribution is the same for 
#' the field and nursery. Lets create functions that will 
#' calculate the likelihood and the negative log-likelihood. 
  Like.s<-function(lambda, dat){
    prod(dpois(dat$slugs, lambda))
  }
  
  minus.logL.s<-function(lambda, dat){
     -sum(dpois(dat$slugs, lambda, log=TRUE))
  }

  
#' Now, lets look at different values of lambda ranging from 0.1 to 5
 lambda.test<-seq(0.1, 5, length=1000)
 nvals<-length(lambda.test) # number of lambdas to explore
 L.vals<-matrix(NA, nvals, 1)
 for(i in 1:nvals){
   L.vals[i]<-Like.s(lambda.test[i], dat = slugs)
 }
 
#' Instead of using a loop, we can use sapply to calculate the 
#' likelihood for several different values of lambda.
 L.vals<-sapply(1:nvals, FUN=function(x){Like.s(lambda.test[x], dat = slugs)}) 
 minus.logL.vals<-sapply(1:nvals, FUN=function(x){minus.logL.s(lambda.test[x], slugs)}) 
  
#' ## Finding the MLE numerically 
#'
#'  Use which.max and which.min to pull of the value of lambda 
#'  that maximizes L(lambda|Y) and  minimizes -logL(lambda | Y), 
#'  respectively  
  lambda.test[which.max(L.vals)] # Gives us lambda that gives the maximum of L.vals
  lambda.test[which.min(minus.logL.vals)] # Gives us lambda that minimizes -LogL
  
#' Now, lets plot the Likelihood and log-likelihood as a function of lambda
#+fig.height=8, fig.width=6  
par(mfcol=c(2,1), bty="L")
  plot(lambda.test, L.vals, xlab=expression(lambda), ylab="Likelihood",type="l")
  abline(v=lambda.test[which.max(L.vals)])
  
  plot(lambda.test, minus.logL.vals, xlab=expression(lambda), ylab="Log Likelihood", type="l")
  abline(v=lambda.test[which.min(minus.logL.vals)])
  
#' Notice: the sample mean is the maximum likelihood estimator in this case!
  mean(slugs$slugs)
  
  
#' ## Now, using optim to get a more precise answer
#' 
#' Optim can be used to find the parameters that minimize any function.  It takes
#' 
#' - a vector of "starting" parameters (first guesses) as its first argument. 
#' - a function to minimize as its second argument (the first argument of this function must correspond to the vector
#'   of parameters that we want to estimate)
#' - you can send other "stuff" (like the data, here)
#' - You can choose from many different methods (here, we will use BFGS)
#' - We will specify that we want it to also return the Hessian (see Textbook)
#'
#' Note: we may get some warnings here about NaN's being produced.  The parameter lambda of the poisson 
#' distribution has to be positive.  If the method "tries" out negative values when searching for the optim,
#' it will result in NaN's and thus, the warnings.
  mle.fit<-optim(par=2, fn=minus.logL.s,  method="BFGS", hessian=T, dat=slugs)  
  mle.fit 
  
#' Notes on the values returned:
#' 
#' - convergence = 0 suggests we successfully found a minimum
#' - value = the -logL value at the minimum (so -value is the logL at the maximum)
  
#' We can get the variance by calculating the inverse of the Hessian using "solve"  
  (SE<-sqrt(solve(mle.fit$hessian))) #SE of lambda^

#' 95% confidence interval based on asymptotic normality of MLEs
  mle.fit$par + c(-1.96, 1.96)*SE
  
#' What about using glm?  We can get the same MLE and its SE (we will learn about glms when we get to generalized
#' linear models).
  glm.fit<-glm(slugs~1, data=slugs, family=poisson())  
  summary(glm.fit)
  exp(glm.fit$coef)
  predict(glm.fit, se=T, type="resp")


#' ## Separate lambdas for the two fields
#' 
#' Create log-likelihood function
  minus.logL.2<-function(lambdas, dat){
    lambda.rookery<-lambdas[1]
    lambda.nursery<-lambdas[2]
    -sum(dpois(dat$slugs, lambda.rookery*I(dat$field=="Rookery") + 
                 lambda.nursery*I(dat$field == "Nursery"),log=TRUE))
    
  }

  minus.logL.2<-function(params, dat){
    beta0<-params[1]
    beta1<-params[2]
    log.lambda<-beta0 + beta1*I(dat$field=="Rookery")
    lambdas<-exp(log.lambda)
    -sum(dpois(dat$slugs, lambdas,log=TRUE))
    
  }
#' Now, estimate the lambdas using optim
  mle.fit.2<-optim(par=c(2,2), fn=minus.logL.2, dat=slugs, method="BFGS", hessian=T)  
  mle.fit.2  
  
#' Likelihood Ratio Test: 2*(LogL[full]- logL[reduced]), which we compare to a chisq with 1 df
  (LR.stat<-2*(-mle.fit.2$value  +mle.fit$value))
  1-pchisq(LR.stat, df=1) 


#' Compare to results using glm
  full.glm<-glm(slugs~field, data=slugs, family=poisson())  
  reduced.glm<-glm(slugs~1, data=slugs, family=poisson())
  anova(full.glm, reduced.glm, test="Chisq")

  

  
#' ## Profile-likelihood Intervals when estimating a single lambda
#' 
#' Recall, our estimate and CI from before: 
#' 
mle.fit$par + c(-1.96, 1.96)*SE

#' Lets find a 95% interval by inverting the LR test for:
#' 
#'  - HO: $\lambda=\lambda_0$
#'  - HA: $\lambda \ne \lambda_0$
#'  
#' Any values of $\lambda_0$ for which 2$[logL[\hat{\lambda}] - logL[\lambda_0]] < \chi^2_1(0.95)$  would NOT be 
#' rejected. 
  qchisq(0.95, 1)/2
  ind<-I(-mle.fit$value +minus.logL.vals < 1.92) # finds all values of lambda that meet this condition
  (CI.95<-range(lambda.test[ind])) # Gives the 95% CI
  
#' Now, lets plot the results!  
 plot(lambda.test, -minus.logL.vals, type='l', xlab=expression(lambda), ylab='log-likelihood', xlim=c(1.2,2.3), ylim=c(-183,-176 ))
 crit.val<--mle.fit$value-1.96
 abline(h=crit.val, col='seagreen', lty=2)
 segments(CI.95[1],-184,CI.95[1],crit.val,lty=2,col=2,lwd=2)
 segments(CI.95[2],-184,CI.95[2],crit.val,lty=2,col=2,lwd=2)
 arrows(CI.95[1],-182.2,CI.95[2],-182.25,length=.1,angle=45,code=3)
 text(1.8,-182.1,expression(paste('confidence interval for ', lambda)), pos=3,cex=.9, col=4)
 arrows(1.775,-178.759,1.775,-176.8383,length=.1,angle=45,code=3)
 text(2.05,-177.8,'  region satisfying\nlog-likelihood inequality',cex=.9, col=4)    


  
#' ## Footer
sessionInfo() 
   