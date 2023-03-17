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
library(modelsummary) # for creating tables with regression output
library(performance) # for checking model assumptions
library(ggthemes) # for colorblind pallete
library(car) # for residual plots
library(dplyr) # for data wrangling
library(MASS) # for glm.nb
library(DHARMa) # for residual plots
library(Data4Ecologists)
data(slugs)
data("longnosedace")

#' ## Models for slugs data

#' Standard Poisson glm (log link)
fit.pois<-glm(slugs ~ field, data = slugs, family = poisson())
summary(fit.pois)

#' Identity link
fit.pois<-glm(slugs ~ field, data = slugs, family = poisson(link = "identity"))
summary(fit.pois)

#' Log link but with means coding
fit.pois<-glm(slugs ~ field - 1, data = slugs, family = poisson())
summary(fit.pois)


#' ## Longnose Dace Models
#' Linear model
lmdace<-lm(longnosedace~acreage+do2+maxdepth+no3+so4+temp, data=longnosedace)
check_model(lmdace, check = c("linearity", "homogeneity", "qq", "normality"))


#' Standard Poisson glm
glmPdace <- glm(longnosedace ~ acreage + do2 + maxdepth + no3 + so4 + temp, 
              data = longnosedace, family = poisson())
summary(glmPdace)
check_model(glmPdace)

#' ## Interpretation
#' 
#' - Parameters (similar to linear regression) but effects on the 
#' log mean
#' - Rate Ratios = exp(beta) = multiplicative effect on the counts
#' 
#+tabdace, message=FALSE, echo=FALSE 
modelsummary(list("Linear regression" = lmdace, "Poisson regression" = glmPdace),
             estimate  = "{estimate} ({std.error})", statistic = NULL, coef_omit = "SD", gof_omit = ".*",
             title="Parameter estimates (SE) from fitted linear and Poisson regression models.")


#+tabdace2, warning=FALSE, echo=FALSE 
modelsummary(list("Poisson regression" = glmPdace), exponentiate = TRUE, 
             estimate = "{estimate} ({conf.low}, {conf.high})", 
             statistic = NULL, gof_omit = ".*", coef_omit = "SD", 
             title="Incidence Rate Ratios = $exp(\\beta)$ (95 percent confidence intervals).")
#'
#'

#' ## Confidence intervals
#' 
#' Can rely on asymptotic Normality for Maximum Likelihood Estimators 
#' 
# Store coefficients and their standard errors
beta.hats <- coef(glmPdace)
ses <- sqrt(diag(vcov(glmPdace)))
round(cbind(beta.hats-1.96*ses, beta.hats+1.96*ses), 3)

#' Or, can use profile likelihood intervals
confint(glmPdace)

#' When calculating CIs for exp(beta), its best to:
#' 
#' 1. Calculate a CI for beta
#' 2. Exponentiate the confidence limits
round(cbind(exp(beta.hats-1.96*ses), exp(beta.hats+1.96*ses)), 3)


#' ## Model checking/residual diagnostics
#' 
#' We will focus on Pearson Residuals:
#' $$\frac{Y_i - \hat{E}[Y_i | X_i]}{\sqrt{\widehat{Var}[Y_i | X_i]}}$$ 
#' 
#' These should have constant variance if the Poisson assumption is
#' reasonable.
#' 
#' - Plot Pearson residuals versus fitted values and predictors
#+ out.width = "90%", fig.height = 6, fig.width =6 
residualPlots(glmPdace)

#' Can also use check_model
#+ out.width = "60%", fig.height = 9, fig.width =6 
performance::check_model(glmPdace)

#' See also the Dharma package (in book) 


#'
#' ## Goodness of fit test using Pearson residuals
# Number of simulations
  nsims <- 10000

# Set up matrix to hold goodness-of-fit statistics
  gfit.sim <- gfit.obs<-matrix(NA, nsims, 1)
  nobs <- nrow(longnosedace) # number of observations
  
# Extract the estimated coefficients and their asymptotic variance/covariance matrix  
# Use these values to generate nsims new beta's  
  beta.hat <- MASS::mvrnorm(nsims, coef(glmPdace), vcov(glmPdace))

# Design matrix for creating new lambda^'s    
  xmat <- model.matrix(glmPdace)
  for(i in 1:nsims){
    # Generate lambda^'s    
    lambda.hat <- exp(xmat%*%beta.hat[i,])
    
    # Generate new simulated values
    new.y <- rpois(nobs, lambda = lambda.hat)
    
    #' Calculate Pearson Chi-squared statistics
    gfit.sim[i,] <- sum((new.y-lambda.hat)^2/(lambda.hat))
    gfit.obs[i,] <- sum((longnosedace$longnosedace-lambda.hat)^2/lambda.hat)
  }
  # Calculate and return the p-value
  (GOF.pvalue <- mean(gfit.sim > gfit.obs))

#' ## Compare check_model versus our own predictive simulation
  
  g1<-ggplot() 
  nsims=100
  for(i in 1:nsims){
    
    # Generate E[Y|X] and var[Y|X], mu.hat and var.hat below    
    lambda.hat <- exp(xmat%*%beta.hat[i,])
      
    # Generate new simulated values
    new.y<-rpois(nobs, lambda = lambda.hat)
    g1 <- g1 +  geom_line(data=data.frame(new.y), aes(new.y), 
                          stat="density", alpha=0.4, col='red') 
  }  
  g1 + 
    geom_line(data=longnosedace, aes(longnosedace), stat="density", lwd=1.5) +
    xlim(0,300)
  
#' So, the sum of squared residuals were always larger for the observed
#' than simulated data, which provides overwhelming evidence that the
#' poisson model is not a good one...
#'   
#' ## Quasi-Poisson (relax the assumption that Var = Mean)
glmQdace<-glm(longnosedace ~ acreage + do2 + maxdepth + no3 + so4 + temp, 
              data = longnosedace, family = quasipoisson())  
summary(glmQdace)

#'
#' ## Negative Binomial Regression

glmNBdace <- MASS::glm.nb(longnosedace ~ acreage + do2 + maxdepth + no3 + so4 + temp, 
                          data = longnosedace)
summary(glmNBdace)


#' ## Residual plots for negative binomial model
#+ out.width="60%",  fig.width=4, fig.height=4
residualPlot(glmPdace, main = "Poisson model")
residualPlot(glmNBdace, main = "Negative Binomial model")

#+ out.width = "90%", fig.height = 9, fig.width =6 
check_model(glmNBdace) 
 

#' ## Goodness of fit test for negative binomial model
# Number of simulations
  nsims<-10000

# Set up matrix to hold goodness-of-fit statistics
  gfit.sim <- gfit.obs <- matrix(NA, nsims, 1)
  nobs<-nrow(longnosedace) # number of observations
  
# Extract the estimated coefficients and their asymptotic variance/covariance matrix  
# Use these values to generate nsims new beta's and thetas 
  beta.hat <- MASS::mvrnorm(nsims, coef(glmNBdace), vcov(glmNBdace))
  theta.hat <- rnorm(nsims, glmNBdace$theta,  glmNBdace$SE.theta)

# Design matrix for creating new lambda^'s    
  xmat <- model.matrix(glmNBdace)
  for(i in 1:nsims){

    # Generate E[Y|X] and var[Y|X], mu.hat and var.hat below    
    mu.hat <- exp(xmat%*%beta.hat[i,])
    var.hat <- mu.hat + mu.hat^2/theta.hat[i]
      
    # Generate new simulated values
    new.y<-rnbinom(nobs, mu = mu.hat, size = theta.hat[i])
   
    # Calculate goodness-of-fit statistics
    gfit.sim[i,]<-sum((new.y - mu.hat)^2/(var.hat))
    gfit.obs[i, ] <- sum((longnosedace$longnosedace - mu.hat)^2/var.hat)
  }
  # Calculate p-value
  (GOF.pvalue <- mean(gfit.sim > gfit.obs))

#' Can also plot the fit statistics for simulated and real data  
#+out.width="60%", fig.width=8, fig.height=5
  fitstats <- data.frame(GOF.sim=gfit.sim, 
                         GOF.obs = gfit.obs,
                         GOF.sim_GOF.obs = gfit.sim-gfit.obs)
  WVPlots::ShadedDensity(frame = fitstats, 
                         xvar = "GOF.sim_GOF.obs",
                         threshold = 0,
                         title = "Distribution: Chi squared (sim data)- Chi squared (obs data)",
                         tail = "right")+
    annotate("text", x=60, y = 0.005, 
             label="Better fit \n to observed data")+
    annotate("text", x=-100, y = 0.005, 
             label="Better fit to \n simulated data") 
  
#' ## Compare check_model versus our own predictive simulation
  
  g1<-ggplot() 
  nsims=100
  for(i in 1:nsims){
    
    # Generate E[Y|X] and var[Y|X], mu.hat and var.hat below    
    mu.hat <- exp(xmat%*%beta.hat[i,])
    var.hat <- mu.hat + mu.hat^2/theta.hat[i]
    
    # Generate new simulated values
    new.y<-rnbinom(nobs, mu = mu.hat, size = theta.hat[i])
    g1 <- g1 +  geom_line(data=data.frame(new.y), aes(new.y), 
                          stat="density", alpha=0.4, col='red') 
  }  
  g1 + 
    geom_line(data=longnosedace, aes(longnosedace), stat="density", lwd=1.5) +
    xlim(0,300)
  
  
#' ## Model comparisons
AIC(glmPdace, glmNBdace)
  
#' We could calculate these comparisons by hand 
-2*sum(dpois(longnosedace$longnosedace, 
               lambda = predict(glmPdace, type = "resp"), 
               log = TRUE)) + 2 * 7
  
-2*sum(dnbinom(longnosedace$longnosedace, 
                 mu = predict(glmNBdace, type = "resp"), 
                 size = glmNBdace$theta, log = TRUE)) + 2 * 8
  

#' Type III tests 
Anova(glmNBdace)

#' Stepwise selection 
MASS::stepAIC(glmNBdace)

 
