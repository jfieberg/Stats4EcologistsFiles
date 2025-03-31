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
#' ## Logistic regression: Bayesian implementation
    
  
## ---------------------------------------------------------------
lrmod<-function(){
  
  # Priors
  beta.voc ~ dt(0, pow(2.5, -2), 1)
  
  # Likelihood
  for(i in 1:n){
     
  }
}


# Bundle data
voc.scaled <- (exp.m$voc - mean(exp.m$voc)) / (0.5 * sd(exp.m$voc))
yearind <- exp.m$year -2004
jagsdata <- list(  )


# Parameters to estimate
params <- c("alpha", "beta.voc", "p")

 

out.p <- jags.parallel(data = jagsdata, parameters.to.save = params, 
                       model.file = lrmod, n.thin= 2, n.chains = 3, 
                       n.burnin = 1000, n.iter = 5000)
 
#' Compare to glm 
MCMCsummary(out.p, params = c("alpha", "beta.voc"), round = 5)
exp.m$voc.scaled<-voc.scaled
exp.m$fyear <- as.factor(exp.m$year)

mod1b <- glm(exp.m$observed ~ fyear + voc.scaled -1, family = binomial(),
             data = exp.m)
cbind(coef(mod1b), confint(mod1b))
 