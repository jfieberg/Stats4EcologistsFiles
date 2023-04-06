#' ---
#' title: "17-LinearMixedEffectsModels3.R"
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
data("RIKZdatdat")


jags.lme<-function(){   
  
  # Priors for the intercepts   
  for (i in 1:n.groups){        
    b0[i] ~ dnorm(0, tau.b0)   # Random intercepts 
  } 
  
  # Hyper priors that describe the distribution of the site-specific intercepts   
  beta0 ~ dnorm(0, 0.001)      # Mean intercept
  tau.b0 <- 1 / (sigma.b0 * sigma.b0)    
  sigma.b0 ~ dunif(0, 100)     # SD of the random intercepts   
  
  # Priors for fixed effect regression parameters 
  beta1 ~ dnorm(0, 0.001)           # Common slope agec 
  tau.eps <- 1 / ( sigma.eps * sigma.eps)       # Residual precision    
  sigma.eps ~ dunif(0, 100)         # Residual standard deviation   
  
  # Calculate site specific intercepts
  for (i in 1:n.groups){        
    beta0i[i] <- beta0 + b0[i] 
  } 
  
  # Likelihood  
  for (i in 1:nobs) {   
    dbh[i] ~ dnorm(mu[i], tau.eps)     # The random variable   
    mu[i] <- beta0 + b0[site[i]] + beta1*agec[i]  # Expectation   
  } 
}


## -------------------------------------------------------------------------------------
jags.lme.alt<-function(){   
  
  # Priors for the intercepts   
  for (i in 1:n.groups){        
    alpha[i] ~ dnorm(beta0, tau.b0)   # Random intercepts 
  } 
  
  # Hyper priors that describe the distribution of the site-specific intercepts   
  beta0 ~ dnorm(0, 0.001)      # Mean intercept
  tau.b0 <- 1 / (sigma.b0 * sigma.b0)    
  sigma.b0 ~ dunif(0, 100)     # SD of the random intercepts   
  
  # Priors for fixed effect regression parameters 
  beta1 ~ dnorm(0, 0.001)           # Common slope agec 
  tau.eps <- 1 / ( sigma.eps * sigma.eps)       # Residual precision    
  sigma.eps ~ dunif(0, 100)         # Residual standard deviation   
  
  # Likelihood  
  for (i in 1:nobs) {   
    dbh[i] ~ dnorm(mu[i], tau.eps)     # The random variable   
    mu[i] <- alpha[site[i]] + beta1*agec[i]  # Expectation   
  } 
}     


## ----jagspineri, cache=TRUE-----------------------------------------------------------
# Bundle data   
jags.data <- list(dbh = pines$dbh, agec = as.numeric(pines$agec), 
                  site = as.numeric(as.factor(pines$site)),
                  n.groups=length(unique(pines$site)),  
                  nobs = nrow(pines))   


# Parameters to estimate    
parameters <- c("beta0i", "beta0", "beta1", "sigma.b0", "sigma.eps")    


# Start Gibbs sampling  
out.ri <- jags.parallel(data=jags.data,  
                        parameters.to.save=parameters, 
                        model=jags.lme, 
                        n.thin=1, n.chains=3, n.burnin=1000, n.iter=10000)   


## -------------------------------------------------------------------------------------
MCMCsummary(out.ri, params=c("beta0", "beta1", "sigma.b0", "sigma.eps"), round=3)
summary(lmer.ri)


## -------------------------------------------------------------------------------------
MCMCsummary(out.ri, params="beta0i", round=3)


## -------------------------------------------------------------------------------------
denplot(out.ri, parms=c("sigma.b0", "sigma.eps"))


## ----jagsrc1, cache=TRUE--------------------------------------------------------------
jags.lme.rc<-function(){
  
  # Priors for the intercepts
  for (i in 1:n.groups){    
    # To allow for correlation between alpha[i] and beta1[i], we need to model
    # their joint (multivariate) distribution 
    beta0i[i] <- B[i,1] # Random intercepts
    beta1i[i] <- B[i,2] # Random slopes for age
    B[i,1:2]~ dmnorm(B.hat[i,], Tau.B[,]) # distribution of the vector (beta0[i], beta1[i])
    B.hat[i,1]<-beta0 # mean for intercepts
    B.hat[i,2]<-beta1 # mean for slopes
  }
  
  # Hyperpriors for intercepts and slopes
  beta0 ~ dnorm(0, 0.001)       
  beta1 ~ dnorm(0, 0.001)
  
  # Hyperpriors for Sigma =  var/cov matrix of the slope/intercept parameters
  sigma.b0 ~ dunif(0,100) # sd intercepts
  sigma.b1 ~ dunif(0,100) # sd of slopes
  cor.b0.b1 ~ dunif(-1,1) # correlation among intercepts and slopes
  Sigma.B[1,1]<-pow(sigma.b0,2)
  Sigma.B[2,2]<-pow(sigma.b1,2)
  Sigma.B[1,2]<-cor.b0.b1*sigma.b0*sigma.b1
  Sigma.B[2,1]<-Sigma.B[1,2]
  
  # Tau = inverse of Sigma (analogous to precision for univariate normal distribution)
  Tau.B[1:2,1:2]<-inverse(Sigma.B[,])
  
  # Prior for within-site errors 
  tau.eps <- 1 / ( sigma.eps * sigma.eps)       # Residual precision
  sigma.eps ~ dunif(0, 100)         # Residual standard deviation
  
  # Likelihood
  for (i in 1:nobs) {
    dbh[i] ~ dnorm(mu[i], tau.eps)  
    mu[i] <- beta0i[site[i]] + beta1i[site[i]]*agec[i]  
  }
} 


# Parameters to estimate
parameters <- c("beta0", "beta1", "betaoi", "beta1i", "sigma.b0", 
                "sigma.b1", "sigma.eps", "cor.b0.b1")


# Start Gibbs sampling
out.rc <- jags.parallel(data = jags.data,
                        parameters.to.save=parameters, 
                        model=jags.lme.rc, 
                        n.thin=1, n.chains=3, n.burnin=100, n.iter=10000) 


## -------------------------------------------------------------------------------------
MCMCsummary(out.rc, params=c("beta0", "beta1", "sigma.b0", 
                             "sigma.b1", "cor.b0.b1", "sigma.eps"), round=3)
summary(lmer.rc)


## -------------------------------------------------------------------------------------
lmer.rc.ind <- lmer(dbh ~ agec + (1 | site) + (0 + agec | site), data=pines)
summary(lmer.rc.ind)

