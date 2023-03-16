#' ---
#' title: "Count based regression models"
#' author: "John Fieberg"
#' date: "2023-03-15"
#' output: 
#'    html_document:
#'      toc: true
#'      toc_depth: 2
#'      toc_float: true
#' ---

#' ## Preamble
#' 
#' Load packages, etc
#+ warning=FALSE,message=FALSE 
library(ggplot2)
theme_set(theme_bw())
knitr::opts_chunk$set(fig.align='center', out.width = "55%", fig.height = 5, fig.width =6)
set.seed(12111) 
library(kableExtra) # for tables
options(kableExtra.html.bsTable = T)
library(ggthemes) # for colorblind pallete
 library(dplyr) # for data wrangling
library(R2jags) # for fitting models using JAGS
library(MCMCvis) # for summarizing JAGS output
library(mcmcplots) # for denplots and traplots
library(Data4Ecologists)
data(slugs)
data("longnosedace")


#' ## Fit Poisson model using JAGS  
pois.fit<-function(){
  
  # Priors for regression parameters
  for(i in 1:np){
    beta[i] ~ dnorm(0,0.001)
  }  
  for(i in 1:n){
    log(lambda[i]) <- inprod(beta[], x[i,])
    longnosedace[i] ~ dpois(lambda[i])
    
    # Fit assessments
    Presi[i] <- (longnosedace[i] - lambda[i]) / sqrt(lambda[i]) # Pearson residuals
    dace.new[i] ~ dpois(lambda[i])    # Replicate data set
    Presi.new[i] <- (dace.new[i] - lambda[i]) / sqrt(lambda[i]) # Pearson residuals
    D[i] <- pow(Presi[i], 2)
    D.new[i] <- pow(Presi.new[i], 2)
  }
  
  # Add up discrepancy measures
  fit <- sum(D[])
  fit.new <- sum(D.new[])
} 


# Bundle data, here we will use the design matrix from R 
xmat <- model.matrix(~acreage + do2 + maxdepth + no3 + so4 + temp, 
                     data = longnosedace)

jagsdata <- list(x = xmat, np = ncol(xmat),
                 n = nrow(longnosedace), 
                 longnosedace = longnosedace$longnosedace)


# Parameters to estimate
params <- c("beta", "lambda", "Presi", "fit", "fit.new", "dace.new")


# Start Gibbs sampler
out.pois <- jags.parallel(data = jagsdata, parameters.to.save = params, 
                          model.file = pois.fit, n.thin = 2, n.chains = 3, n.burnin = 5000, 
                          n.iter = 20000)


#' Look at results for the regression coefficients and also make sure
#' that out sampler has converged to the same place for all 3 chains
MCMCsummary(out.pois, params="beta", round=3)
denplot(out.pois, parms="beta")

#' Let's plot Pearson residuals versus fitted values
#' for the Poisson regression model
#+  out.width="55%"
bresids <- MCMCpstr(out.pois, params = "Presi", func = mean)$Presi
bfitted <- MCMCpstr(out.pois, params = "lambda", func = mean)$lambda
jagsPfit <- data.frame(resid = bresids, fitted = bfitted)
ggplot(jagsPfit, aes(fitted, resid)) + geom_point() + 
  geom_hline(yintercept = 0) + geom_smooth() +
  xlab("Fitted value") + ylab("Pearson Residual")


#' GOF test 
fitstats <- MCMCpstr(out.pois, params = c("fit", "fit.new"), type = "chains") 
T.extreme <- fitstats$fit.new >= fitstats$fit
(p.val <- mean(T.extreme))


#' ## Fit negative binomial model using JAGS  
nb.fit<-function(){
  
  # Priors for regression parameters
  for(i in 1:np){
    beta[i] ~ dnorm(0,0.001)
  }  
  theta~dunif(0,50)
  for(i in 1:n){
    log(mu[i]) <- inprod(beta[], x[i,])
    p[i] <- theta/(theta+mu[i])
    longnosedace[i] ~ dnegbin(p[i],theta)
    
    # Fit assessments
    Presi[i] <- (longnosedace[i] - mu[i]) / sqrt(mu[i]+mu[i]*mu[i]/theta) # Pearson residuals
    dace.new[i] ~ dnegbin(p[i],theta)    # Replicate data set
    Presi.new[i] <- (dace.new[i] - mu[i]) / sqrt(mu[i]+mu[i]*mu[i]/theta) # Pearson residuals
    D[i] <- pow(Presi[i], 2)
    D.new[i] <- pow(Presi.new[i], 2)
  }
  
  # Add up discrepancy measures
  fit <- sum(D[])
  fit.new <- sum(D.new[])
} 


# Parameters to estimate
params <- c("mu", "beta", "theta","Presi", "fit", "fit.new", "dace.new")


# Start Gibbs sampler
out.nb <- jags.parallel(data = jagsdata, parameters.to.save = params, 
                        model.file = nb.fit, n.thin = 2, n.chains = 3, n.burnin = 5000, 
                        n.iter = 20000)


#' Let's again plot Pearson residuals versus fitted values for 
#' the Negative Binomial regression model 
#+ out.width="55%"
nbresids <- MCMCpstr(out.nb, params = "Presi", func = mean)$Presi
nbfitted <- MCMCpstr(out.nb, params = "mu", func = mean)$mu
jagsnbfit <- data.frame(resid = nbresids, fitted = nbfitted)
ggplot(jagsnbfit, aes(fitted, resid)) + geom_point() + 
  geom_hline(yintercept = 0) + geom_smooth() +
  xlab("Fitted value") + ylab("Pearson Residual")


# GOF test
fitstats <- MCMCpstr(out.nb, params = c("fit", "fit.new"), type = "chains") 
T.extreme <- fitstats$fit.new >= fitstats$fit
(p.val <- mean(T.extreme))


#' ## Comparison of frequentist and Bayesian Negative Binomial Models 

#' Frequentist fit
data("longnosedace")
glmNBdace <- MASS::glm.nb(longnosedace ~ acreage + do2 + maxdepth + no3 + so4 + temp, 
                          data = longnosedace)

# Summary of the Bayesian model
Nbsum<-MCMCsummary(out.nb, params=c("beta","theta"), round=3)

# Confidence intervals Frequentist model
ci.nb<-confint(glmNBdace)

# Create data frame with results
est.dat<-data.frame(
  bhat = c(Nbsum[, 1], coef(glmNBdace), glmNBdace$theta),
  lowci = c(Nbsum[, 3], ci.nb[, 1], glmNBdace$theta - 1.96 * glmNBdace$SE.theta),
  upci = c(Nbsum[, 5], ci.nb[, 2], glmNBdace$theta + 1.96 * glmNBdace$SE.theta),
  Method = c(rep("Freq", 8), rep("Bayes", 8)),
  coefname = rep(c(names(coef(glmNBdace)), "theta"), 2))

ggplot(est.dat, aes(coefname, bhat, colour = Method)) + 
  geom_point(size = 5, shape = 18, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymax = upci, ymin = lowci, colour = Method),
                position = position_dodge(width = 0.5)) + ylab("Point Estimate") +
  xlab("Coefficient") +  scale_colour_colorblind()


## ---------------------------------------------------------------------------------------------
#' Fit Poisson-normal model using JAGS  
poisnorm.fit<-function(){
  
  # Priors for regression parameters
  for(i in 1:np){
    beta[i] ~ dnorm(0,0.001)
  }  
  
  # Prior for sigma and tau
  sig ~ dunif(0, 3) 
  tau <- 1/(sig*sig)
  
  
  for(i in 1:n){
    # Likelihood
    log(lambda[i]) <- inprod(beta[], x[i,]) + eps[i]
    longnosedace[i] ~ dpois(lambda[i])
    eps[i] ~ dnorm(0, tau)
    
    # Calculate mean and variance
    mu[i] <-  exp(inprod(beta[], x[i,]) + sig^2/2)
    vary[i] <- mu[i] + (exp(sig^2)-1)*mu[i]^2
    
    # Fit assessments
    Presi[i] <- (longnosedace[i] - mu[i]) / sqrt(vary[i]) # Pearson residuals
    dace.new[i] ~ dpois(lambda[i])    # Replicate data set
    Presi.new[i] <- (dace.new[i] - mu[i]) / sqrt(vary[i]) # Pearson residuals
    D[i] <- pow(Presi[i], 2)
    D.new[i] <- pow(Presi.new[i], 2)
  }
  
  # Add up discrepancy measures
  fit <- sum(D[])
  fit.new <- sum(D.new[])
} 

# Parameters to estimate
params <- c("beta", "mu", "Presi", "fit", "fit.new")


# Start Gibbs sampler
out.poisnorm <- jags.parallel(data = jagsdata, parameters.to.save = params, 
                              model.file = poisnorm.fit, n.thin = 2, n.chains = 3, n.burnin = 5000, 
                              n.iter = 20000)


#' Pearson residuals versus fitted values for the Poisson-Normal regression model
#+ out.width="55%" 
pnresids <- MCMCpstr(out.poisnorm, params = "Presi", func = mean)$Presi
pnfitted <- MCMCpstr(out.poisnorm, params = "mu", func = mean)$mu
jagspnfit <- data.frame(resid = pnresids, fitted = pnfitted)
ggplot(jagspnfit, aes(fitted, resid)) + geom_point() + 
  geom_hline(yintercept = 0) + geom_smooth() +
  xlab("Fitted value") + ylab("Pearson Residual")


#' GOF 
fitstats <- MCMCpstr(out.poisnorm, params = c("fit", "fit.new"), type = "chains") 
T.extreme <- fitstats$fit.new >= fitstats$fit
(p.val <- mean(T.extreme))


#' Comparisons using DIC 
out.pois$BUGSoutput$DIC
out.nb$BUGSoutput$DIC
out.poisnorm$BUGSoutput$DIC

