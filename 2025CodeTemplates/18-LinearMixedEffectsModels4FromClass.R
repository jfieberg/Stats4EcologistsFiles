#' ---
#' title: "17-LinearMixedEffectsModels4.R"
#' author: "John Fieberg"
#' output: 
#'    html_document:
#'      toc: true
#'      toc_depth: 2
#'      toc_float: true
#' ---

#' ## Load libraries
#+ warning=FALSE,message=FALSE
library(ggplot2);theme_set(theme_bw())
knitr::opts_chunk$set(fig.align='center', 
                      out.width = "55%", 
                      fig.height = 5, 
                      fig.width =6,
                      error = FALSE)
library(tidyverse) # for data wrangling
library(gridExtra) # for multi-panel plots
library(lme4)  # for fitting random-effects models
library(nlme) # for fitting random-effects models
library(glmmTMB) # for fitting random-effects models
library(sjPlot) # for visualizing fitted models
library(modelsummary) # for creating output tables
library(kableExtra) # for tables 
options(kableExtra.html.bsTable = T)
library(ggplot2)# for plotting
library(performance) # for model diagnostics 
library(R2jags) # to call JAGS from R
library(MCMCvis) # summarize JAGS output
library(mcmcplots)# summarize JAGS output

#' Load data sets
#+warning=FALSE, message=FALSE 
library(Data4Ecologists) # for data
data("RIKZdat")

#' Specify the model similar to lmer (in terms of betas and b's)
#' How many priors do we need?

jags.lme<-function(){   
  
  # Priors for the b's
  for (i in 1:9){        
    b0[i] ~ dnorm(0, tau.b0)
  } 
  
  # Priors for the betas
  beta0 ~ dnorm(0, 0.001)
  beta1 ~ dnorm(0, 0.001)
  beta2 ~ dnorm(0, 0.001)
  
  # Priors for the variance parameters
  sigma ~ dunif(0, 10)
  sigma.b0 ~ dunif(0, 10)
  tau<-1/(sigma*sigma)
  tau.b0 <- 1/(sigma.b0*sigma.b0)
  
  # Calculate site specific intercepts
  for (i in 1:9){        
     gamma[i]<- beta0 + b0[i]
  } 
  
  # Likelihood  
  for (i in 1:nobs) {   
    Richness[i] ~ dnorm(mu[i], tau)
    mu[i] <- beta0 + b0[Beach[i]] + beta1*NAPc[i] 
  } 
}

# Bundle data   
jags.data <- list(nobs = nrow(RIKZdat), 
                  NAPc = as.numeric(scale(RIKZdat$NAP)), 
                  Beach = RIKZdat$Beach, 
                  Richness = RIKZdat$Richness)

# Parameters to estimate    
parameters <- c("beta0", "beta1", "beta2", "sigma", "sigma.b0" )    


# Start Gibbs sampling  
out.ri <- jags.parallel(data=jags.data,  
                        parameters.to.save=parameters, 
                        model=jags.lme, 
                        n.thin=1, n.chains=3, n.burnin=1000, n.iter=10000)   


## -------------------------------------------------------------------------------------

#' Specify the model in terms of $\beta_{0i} = \beta_0 + b_{0i}$

jags.lme<-function(){   
  
  # Priors for the b's
  for (i in 1:9){        
    gamma[i] ~ dnorm(beta0, tau.b0)
  } 
  
  # Priors for the betas
  beta0 ~ dnorm(0, 0.001)
  beta1 ~ dnorm(0, 0.001)
  beta2 ~ dnorm(0, 0.001)
  
  # Priors for the variance parameters
  sigma ~ dunif(0, 10)
  sigma.b0 ~ dunif(0, 10)
  tau<-1/(sigma*sigma)
  tau.b0 <- 1/(sigma.b0*sigma.b0)
  
  # Likelihood  
  for (i in 1:nobs) {   
    Richness[i] ~ dnorm(mu[i], tau)
    mu[i] <- gamma[Beach[i]] + beta1*NAPc[i] 
  } 
}

# Bundle data   
jags.data <- list(nobs = nrow(RIKZdat), 
                  NAPc = as.numeric(scale(RIKZdat$NAP)), 
                  Beach = RIKZdat$Beach, 
                  Richness = RIKZdat$Richness)

# Parameters to estimate    
parameters <- c("beta0", "beta1", "beta2", "sigma", "sigma.b0", "gamma" )    


# Start Gibbs sampling  
out.ri <- jags.parallel(data=jags.data,  
                        parameters.to.save=parameters, 
                        model=jags.lme, 
                        n.thin=1, n.chains=3, n.burnin=1000, n.iter=10000)   
