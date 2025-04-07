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
  for (i in 1:n.groups){        
   
  } 
  
  # Priors for the betas
 
  
  # Priors for the variance parameters
  
  
  # Calculate site specific intercepts
  for (i in 1:n.groups){        
     
  } 
  
  # Likelihood  
  for (i in 1:nobs) {   
  
    
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


## -------------------------------------------------------------------------------------

#' Specify the model in terms of $\beta_{0i} = \beta_0 + b_{0i}$
jags.lme.alt<-function(){   
  
  # Priors for the intercepts   
  for (i in 1:n.groups){        

  } 
  
  # Hyper priors that describe the distribution 
  # of the site-specific intercepts (their mean and variance)  
   
  
  # Priors for other fixed effect regression parameters (exposure)

  
  # Priors for within beach residuals 
  
  
  # Likelihood  
  for (i in 1:nobs) {   

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


 
