#' ---
#' title: "Logistic regression models"
#' author: "John Fieberg"
#' output: 
#'    html_document:
#'      toc: true
#'      toc_depth: 2
#'      toc_float: true
#' ---
#' 
## ----echo=FALSE, warning=FALSE,message=FALSE--------------------
set.seed(12111)
library(ggplot2) # for plotting
library(MASS) # for generating multivariate random normal variables
library(R2jags) # for fitting models using JAGS
library(MCMCvis) # for visualizing MCMC output
library(mcmcplots) # for denplot and traplot for inspecting MCMC output
library(ggthemes) # for colorblind pallete
knitr::opts_chunk$set(fig.align = 'center',
                      out.width = "55%",
                      fig.height = 5,
                      fig.width = 6)

#' Read in a moose data set  
library(SightabilityModel)
data(exp.m)
str(exp.m)

#'  
#' ## Logistic regression: Bayesian implementations
    
#' We need to think a bit abour prior distributions and their implications
#' for probabilites after backtransforming from the logit scale. 
#+fig.align="center",   fig.height=4, fig.width=8, out.width = "75%"
par(mfrow=c(1,2))
beta0<-rnorm(10000,0, 1/ sqrt(0.001))
p<-plogis(beta0)
hist(beta0, xlab=expression(beta[0]), 
     main=expression(paste("Prior distribution for ", beta[0])), 
     col="gray")
hist(p, xlab="p", main="Prior distribution for p", col="gray")


#' 
#' N(0, tau = 1/3) is probably a better option  
#'  
#+ fig.align="center", fig.height=4, fig.width=8, out.width = "75%"
par(mfrow=c(1,2))
beta0<-rnorm(10000, 0, sqrt(3))
p<-plogis(beta0)
hist(beta0, xlab=expression(beta[0]), 
     main=expression(paste("Prior distribution for ", beta[0])), 
     col="gray")
hist(p, xlab="p", main="Prior distribution for p", col="gray")


#' 
#' Gelman's recommendations: 
#' 
#' 1. Scaling continuous predictors so they have mean 0 and sd = 0.5.
#' 2. Using a Cauchy prior, with precision parameter = 0.16
#' 
#' This distribution is equivalent to a Student-t distribution with 1 degree 
#' of freedom and can be specified in JAGS as: `dt(0, pow(2.5,-2), 1)`. 
#'  
#'  
#+ fig.align="center", fig.keep="last", out.width="50%",  message=FALSE, warning=FALSE
library(LaplacesDemon)
curve(dnorm(x, mean=0, sd=2.5), from=-8, to=8, ylab="Density")
curve(dstp(x, mu=0, tau=1/2.5^2, n=1), from=-8, to=8, add=TRUE, lty=2)
legend(2, 0.15, c("Normal", "t"), lty=c(1,2), bty="n")


#' 
#' ### Fitting the logistic regression model to moose data
#' 
#' Below, we use Gelman's suggested prior 
#' 
## ---------------------------------------------------------------
lrmod<-function(){
  
  # Priors
  alpha ~ dt(0, pow(2.5, -2), 1)
  beta ~ dt(0, pow(2.5, -2), 1)
  
  # Likelihood
  for(i in 1:n){
    logit(p[i]) <-  
    observed[i]  
    
    # GOF test
    presi[i] <- 
      
    obs.new[i]  
    presi.new[i] <- (obs.new[i] - p[i]) / sqrt(p[i] * (1 - p[i]))
  
    D[i] <- pow(presi[i], 2)
    D.new[i] <- pow(presi.new[i], 2)
  }
  fit <- sum(D[])
  fit.new <- sum(D.new[])
}


# Bundle data
voc.scaled <- (exp.m$voc - mean(exp.m$voc)) / (0.5 * sd(exp.m$voc))
jagsdata <- list(observed = exp.m$observed, voc = voc.scaled, n = nrow(exp.m))


# Parameters to estimate
params <- c("alpha", "beta", "p", "presi", "fit", "fit.new")

# MCMC settings
nc <- 3
ni <- 3000
nb <- 1000
nt <- 2

out.p <- jags.parallel(data = jagsdata, parameters.to.save = params, 
                       model.file = lrmod, n.thin= 2, n.chains = 3, 
                       n.burnin = 1000, n.iter = 3000)

#Goodness-of-fit test
fitstats <- MCMCpstr(out.p, params = c("fit", "fit.new"), type = "chains") 
T.extreme <- fitstats$fit.new >= fitstats$fit
(p.val <- mean(T.extreme))

#' Compare to glm 
MCMCsummary(out.p, params = c("alpha", "beta"), round = 3)
exp.m$voc.scaled<-voc.scaled

mod1b <- glm(exp.m$observed ~ voc.scaled, family = binomial(),
             data = exp.m)
cbind(coef(mod1b), confint(mod1b))

#' 
#' 
#' What if we had naively used extremely vague priors for our regression parameters?  
#' We can compare the results, below, using the `MCMCplot` function 
#' in the `MCMCvis` package:
#' 
## ---------------------------------------------------------------
lrmodv<-function(){
  
  alpha~dnorm(0, 0.0001)
  beta~dnorm(0, 0.0001)
  
  # Likelihood
  for(i in 1:n){
    logitp[i]<-alpha+beta*voc[i]
    p[i]<-exp(logitp[i])/(1+exp(logitp[i]))
    observed[i]~dbin(p[i],1)
  }
}

params <- c("alpha", "beta", "p")
out.p.vague <- jags.parallel(data = jagsdata, parameters.to.save = params, 
                             model.file = lrmodv, n.thin= 2, n.chains = 3, 
                             n.burnin = 1000, n.iter = 3000)

#' Compare the two posterior distributions (using different priors)
MCMCplot(object = out.p, object2=out.p.vague, params=c("alpha", "beta"), 
         offset=0.1, main='Posterior Distributions')

 