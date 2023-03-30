#' ---
#' title: "Zero-inflation "
#' author: "John Fieberg"
#' output: 
#'    html_document:
#'      toc: true
#'      toc_depth: 2
#'      toc_float: true
#' ---

#' Load libraries
## ----warning=FALSE, message=FALSE--------------------------
set.seed(1288)
library(patchwork) # for creating multi-panel plots
library(ggplot2) # for plotting
library(R2jags) # for interfacing with JAGS
library(MCMCvis) # for summarizing MCMC output
library(mcmcplots) # for visualizing MCMC output
library(sjPlot) # for creating tables
library(modelsummary) # for creating tables 
library(dplyr) # for data wrangling
library(knitr) # for creating reproducible reports
library(MASS) # for fitting negative binomial regression models
library(pscl) # for fitting zero-inflated models


## ----echo=FALSE, warning=FALSE,message=FALSE---------------
theme_set(theme_bw())
knitr::opts_chunk$set(
  fig.align = 'center',
  out.width = "55%",
  fig.height = 5,
  fig.width = 6
)


#' Load the data
library(Data4Ecologists)
data(slugs) 


#' ## Goodness of fit test using number of zeros

#' Fitted model
  fit.pois<-glm(slugs~field, data=slugs, family=poisson())

# Number of simulations
  nsims<-10000

# Set up matrix to hold goodness-of-fit statistics
  zeros.sim<-matrix(NA, nsims, 1)
  nobs<-nrow(slugs) # number of observations
  
# Extract the estimated coefficients and their asymptotic variance/covariance matrix  
# Use these values to generate nsims new beta's  
  beta.hat<-MASS::mvrnorm(nsims, coef(fit.pois), vcov(fit.pois))

# Design matrix for creating new lambda^'s    
  xmat<-model.matrix(fit.pois)
  for(i in 1:nsims){
    # Generate lambda^'s    
    lambda.hat<-exp(xmat%*%beta.hat[i,])
    # Generate new simulated values
    new.y<-rpois(nobs, lambda = lambda.hat)
    # Count the number of zeros
    zeros.sim[i]<-sum(new.y==0)
  }


## ----fig.height = 5, fig.width = 6, out.width = "55%"
# T for our data
(T<-sum(slugs$slugs==0))
hist(zeros.sim, xlab="# of zeros", ylab="Density", main="", col="gray")
abline(v=T)

#pvalue
sum(zeros.sim>=T)/nsims


#' ## Fishing example
data(fish)


#' ## Zero-inflated Poisson 
m.zip <- zeroinfl(count ~ child + camper +persons | child, 
                  data = fish, dist="poisson")
summary(m.zip)


#' Calculate the probability of a non-inflated zero (get to go fishing)
#' when you have 3 children 
(     )

#' Calculate the multiplicative effect on the mean number of fish 
#' caught for every 1 unit increase in the number of children in your party. 
( ) 



#' ## Zero-inflated negative binomial
m.zinb <- zeroinfl(count ~ child + camper +persons | child,
               data = fish, dist = "negbin")
summary(m.zinb)


#' Non-zero inflated models 
m.pois<- glm(count ~ child + camper + persons, family = poisson, data = fish)
m.nb <- glm.nb(count ~ child + camper + persons, data = fish)
AIC(m.pois, m.nb, m.zip, m.zinb)

 
#' Pearson residual plots
#+fig.align="center", out.width="80%", fig.height=4, fig.width=8
fish$predict.m.zinb<-predict(m.zinb, type="response")
fish$presiduals.m.zinb<-resid(m.zinb, type="pearson")
p1<-ggplot(fish, aes(predict.m.zinb, count))+geom_point()+
    geom_jitter()+geom_abline(intercept=0,slope=1)+
  xlab("Predicted")+ylab("Observed")

p2<-ggplot(fish, aes(predict.m.zinb, presiduals.m.zinb))+geom_point()+
    geom_jitter()+
    geom_abline(intercept=0,slope=0)+
  xlab("Predicted")+ylab("Pearson Residual")
p1+p2


#' Inspect the large outlier (catch) and large outlier (Pearson resid) 
ind1<-which.max(fish$count) # Max count
ind2<-which.max(fish$presiduals.m.zinb) # Max Pearson residual
fish[c(ind1, ind2),] 


#' ## Bayesian implementation 
#' 
#' First, create an indicator variable for "not inflated zero" (i.e., went fishing)
I.fish<-rep(NA, nrow(fish))
I.fish[fish$count>0]<-1 # these are Not inflated zeros

# Bundle data
jagsdata <- list(count=fish$count, child=fish$child, 
                 camper=as.numeric(fish$camper)-1, 
                 persons=fish$persons, n=nrow(fish), 
                 I.fish=I.fish)

#' JAGS zero-inflated Poisson model
znb<-function(){
  
  # Priors for count model
  for(i in    ){
    beta.c[i] 
  }
  
  # Priors for zero-inflation model
  for(i in    ){
    beta.zi[i] 
  } 
  
  # Overdispersion parameter
  theta 
  
  # Likelihood
  for(i in 1:n){
    # Zero-inflation part of the model
    # Let psi = probability went fishing (not inflated zero) 
      
     
    # Count part of the model
    # Let mu = mean given the party went fishing
    # Let mu.eff = marginal (or unconditional) mean number of fish 
    # caught(pooling data from those that fished or did not fish)
    # mu.eff = mu*I.fish   
     
    # Mean and variances of the observations
    # From eq 11.23 of Zuur et al. (p. 277), but using psi = 1-pi    
    Ey[i]<-mu[i]*psi[i]  
    Vary[i]<-psi[i]*(mu[i])*(1+mu[i]*(1-psi[i])) 
    
    # Generate "new" data
    I.fish.new[i]~dbin(psi[i],1)
    mu.eff.new[i]<-mu[i]*I.fish.new[i]
    count.new[i]~dpois(mu.eff.new[i])
    
    # Pearson residuals
    presi[i]<-(count[i]-Ey[i])/sqrt(Vary[i]) # Pearson Resid
    presi.new[i]<-(count.new[i]-Ey[i])/sqrt(Vary[i])
    
    # Discrepancy measures
    D[i]<-pow(presi[i], 2)
    D.new[i]<-pow(presi.new[i],2)
  }
  fit<-sum(D[])
  fit.new<-sum(D.new[])
}



# Parameters to estimate
params <- c("beta.zi", "beta.c", "Ey", "psi", "mu", "presi",
            "presi.new", "fit", "fit.new", "theta", "I.fish")

# MCMC settings
nc <- 3
ni <- 15000
nb <- 4000
nt <- 10

# Start sampler 
out.znb <- jags.parallel(data = jagsdata, parameters.to.save = params, 
                         model.file = znb, n.thin = 10, n.chains = 3, n.burnin = 4000, 
                         n.iter= 15000)
MCMCsummary(out.znb, params = c("beta.zi", "beta.c", "theta"))




#' JAGS zero-inflated negative binomial model
znb<-function(){
  
# Priors for count model
 for(i in 1:4){
   beta.c[i]~dnorm(0,0.001)
 }
 
# Priors for zero-inflation model
 for(i in 1:2){
   beta.zi[i]~dnorm(0,1/3)
 } 
  
# Overdispersion parameter
  theta~dunif(0,50)
    
# Likelihood
  for(i in 1:n){
  # zero-inflation part (logit prob NOT inflated 0, i.e., "went fishing"))  
    logit(psi[i])<-beta.zi[1]+beta.zi[2]*child[i] 
    I.fish[i]~dbern(psi[i]) # NOt zero-inflated (i.e., "went fishing") 
     
  # Count part
    log.mu[i]<-beta.c[1]+beta.c[2]*child[i] + beta.c[3]*camper[i]+beta.c[4]*persons[i]
    mu[i]<-exp(log.mu[i])
    mu.eff[i]<-mu[i]*I.fish[i]
    p[i]<-theta/(theta+mu.eff[i])
    count[i]~dnegbin(p[i], theta)
  
  # Mean and variances of the observations
  # From eq 11.23 of Zuur et al. (p. 277), but using psi = 1-pi    
    Ey[i]<-mu[i]*psi[i]  
    Vary[i]<-psi[i]*(mu[i]+mu[i]^2/theta)+mu[i]^2*((1-psi[i])^2 +(1-psi[i])) 
  
  # Generate "new" data
    I.fish.new[i]~dbin(psi[i],1)
    mu.eff.new[i]<-mu[i]*I.fish.new[i]
    p.new[i]<-theta/(theta+mu.eff.new[i])
    count.new[i]~dnegbin(p.new[i], theta)
  
  # Pearson residuals
    presi[i]<-(count[i]-Ey[i])/sqrt(Vary[i]) # Pearson Resid
    presi.new[i]<-(count.new[i]-Ey[i])/sqrt(Vary[i])

  # Discrepancy measures
    D[i]<-pow(presi[i], 2)
    D.new[i]<-pow(presi.new[i],2)
  }
  fit<-sum(D[])
  fit.new<-sum(D.new[])
}

# Bundle data
jagsdata <- list(count=fish$count, child=fish$child, 
                 camper=as.numeric(fish$camper)-1, 
                 persons=fish$persons,n=nrow(fish), 
                 I.fish=I.fish )
                 
# Parameters to estimate
params <- c("beta.zi", "beta.c", "Ey", "psi", "mu", "presi",
            "presi.new", "fit", "fit.new", "theta", "I.fish")

# MCMC settings
nc <- 3
ni <- 15000
nb <- 4000
nt <- 10

# Start sampler 
out.znb <- jags.parallel(data = jagsdata, parameters.to.save = params, 
                model.file = znb, n.thin = 10, n.chains = 3, n.burnin = 4000, 
                n.iter= 15000)
MCMCsummary(out.znb, params = c("beta.zi", "beta.c", "theta"))


#' Plot output
#+ out.width = "100%" 
  denplot(out.znb, parms=c("beta.zi", "beta.c", "theta"), ask=FALSE)
  traplot(out.znb, parms=c("beta.zi", "beta.c", "theta"), ask=FALSE)


#' Bayesian p-value
fitstats <- MCMCpstr(out.znb, params = c("fit", "fit.new"), type = "chains") 
T.extreme <- fitstats$fit.new >= fitstats$fit
(p.val <- mean(T.extreme))  


#' Infer who went fishing! 
fish$I.fish.hat<-out.znb$BUGSoutput$mean$I.fish
ggplot(fish, aes(x=as.factor(child), y=I.fish.hat))+
  geom_boxplot()+ xlab("Number of Children")+ylab("Went fishing?")+
   geom_jitter(color="black", size=1, alpha=0.9)


## ----------------------------------------------------------
fish %>% group_by(child) %>% dplyr::summarize(meanfish=mean(count))



#' Equivalent zero-inflated Poisson models using glmmTMB
#+ warning=FALSE, message = FALSE
library(glmmTMB)
m.zip.TMB <- glmmTMB(count ~ child + camper + persons, ziformula = ~ child, 
                     data = fish, family  = poisson)
m.zip
m.zip.TMB


#' Equivalent zero-inflated Negative Binomial Models
m.zinb.TMB <- glmmTMB(count ~ child + camper + persons, ziformula = ~ child,
                    data = fish, family = nbinom2)
m.zinb
m.zinb.TMB


