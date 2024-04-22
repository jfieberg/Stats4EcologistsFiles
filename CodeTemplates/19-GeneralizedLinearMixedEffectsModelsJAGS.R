#' ---
#' title: "19-GeneraizedLinearMixedEffectsModelsJAGS.R"
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

#' Only 1 beach has the lowest exposure level:  Recode exposure to have only 2 levels.	
RIKZdat$exposure.c<-"High"	  
RIKZdat$exposure.c[RIKZdat$exposure%in%c(8,10)]<-"Low"	
RIKZdat$NAPc = RIKZdat$NAP-mean(RIKZdat$NAP) #center NAP variable	


glmer.ri <- glmer(Richness ~ NAPc + exposure.c + (1|Beach),
                 family=poisson(), data=RIKZdat)

#' Specify the model similar to lmer (in terms of betas and b's)
#' How many priors do we need?

jags.lme<-function(){   
  
  # Distribution of the b's
  for (i in 1:n.groups){        
   b[i] ~ dnorm(0, taub)
  } 
  
  # Priors for the betas
  beta0 ~ dnorm(0, 0.0001)
  beta1 ~ dnorm(0, 0.0001)
  beta2 ~ dnorm(0, 0.0001)
  
  # Priors for the variance parameters
  sigmab ~ dunif(0, 5)
  taub <-1/(pow(sigmab, 2))
  
  sigeps ~ dunif(0, 5)
  taueps <- 1/(pow(sigeps, 2))
  
  # Calculate site specific intercepts
  for (i in 1:n.groups){        
    Beach.int[i] <-beta0 + b[i]  
  } 
  
  # Likelihood  
  for (i in 1:nobs) {   
    Richness[i] ~ dnorm(mu[i], taueps)
    mu[i] <- beta0 +  b[Beach[i]] + beta1*NAPc[i] + beta2*exposure[i]
  } 
}

# Bundle data   
jags.data <- list(n.groups = length(unique(RIKZdat$Beach)), 
                  nobs = nrow(RIKZdat), 
                  Richness = RIKZdat$Richness, 
                  NAPc = RIKZdat$NAPc,
                  exposure = ifelse(RIKZdat$exposure.c=="High", 0, 1),
                  Beach = as.numeric(RIKZdat$Beach))

# Parameters to estimate    
parameters <- c("beta0", "beta1", "beta2", "sigmab", "sigeps", "Beach.int")    


# Start Gibbs sampling  
out.ri <- jags.parallel(data=jags.data,  
                        parameters.to.save=parameters, 
                        model=jags.lme, 
                        n.thin=1, n.chains=3, n.burnin=1000, n.iter=10000)   


## -------------------------------------------------------------------------------------

#' Specify the model in terms of $\beta_{0i} = \beta_0 + b_{0i}$
jags.lme.alt<-function(){   
  
  # Beach-level the intercepts   
  for (i in 1:n.groups){        
     Beach.int ~ dnorm(beta0, taub)
  } 
  
  # Priors for the betas
  beta0 ~ dnorm(0, 0.0001)
  beta1 ~ dnorm(0, 0.0001)
  beta2 ~ dnorm(0, 0.0001)
  
  # Priors for the variance parameters
  sigmab ~ dunif(0, 5)
  taub <-1/(pow(sigmab, 2))
  
  sigeps ~ dunif(0, 5)
  taueps <- 1/(pow(sigeps, 2))
  
  # Likelihood  
  for (i in 1:nobs) {   
    Richness[i] ~ dnorm(mu[i], taueps)
    mu[i] <- beta0 +  b[Beach[i]] + beta1*NAPc[i] + beta2*exposure[i]
  } 
}     

# Bundle data   
jags.data <- list( )

# Parameters to estimate    
parameters <- c( )    


# Start Gibbs sampling  
out.ri <- jags.parallel(data=jags.data,  
                        parameters.to.save=parameters, 
                        model=jags.lme, 
                        n.thin=1, n.chains=3, n.burnin=1000, n.iter=10000)   


denplot(out.ri) 
