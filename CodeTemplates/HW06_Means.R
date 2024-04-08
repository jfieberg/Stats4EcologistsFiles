#' ---
#' title: "HW06: Generalized Linear Models"
#' author: ""
#' date: "`r format(Sys.time(), '%d %B, %Y')`"
#' output: pdf_document
#' ---
#' 
#' ## Preamble
#' 
## ----setup, message=FALSE-----------------------------------------------------
# set global chunk options: images will be 7x5 inches
# wrap output at 60 characters 
library(R2jags) # for interfacing with JAGS
library(mcmcplots) # for plotting results
library(ggplot2) # for some later plotting
library(MCMCvis) # for summarizing MCMC results 
library(dplyr) # for data wrangling
library(knitr) # for options when knitting
library(ggplot2) # for plotting
library(ggthemes) # for colorblind theme
library(carData) # to read in data
opts_chunk$set(fig.width=7, fig.height=5, tidy.opts = list(width.cutoff = 60), tidy = TRUE)  

#' 
#' ## 1. Exercise 12.2
#' ##1. Logistic Regression Model #1
#' 
## ----message = FALSE----------------------------------------------------------
# read in data
data("TitanicSurvival")

#' Drop missing values
TitanicSurvival <- TitanicSurvival %>%
  filter(!is.na(age)) %>%
  mutate(age.sc = (age - mean(age))/sd(age),
         survived.ind = ifelse(survived == "yes", 1, 0))


#' Fit model with effects coding
modeffects <- glm(survived.ind ~ sex + age + passengerClass, 
                  family = binomial(), data = TitanicSurvival)
summary(modeffects)

#' 1. Interpret coefficients....
#' 2. Note the number of parameters (5 total)...with 5 coefficients, we can
#' estimate any combination of male/female x passenger class


#' Now, try fitting the model with means coding
modmeans_a <- glm(survived.ind ~ sex + age + passengerClass -1, 
                  family = binomial(), data = TitanicSurvival)
summary(modmeans_a)

#' 1. Interpret coefficients....
#' 2. Note the number of parameters (5 total)

#' Now, version 2, changing the order in which the categorical variables enter
modmeans_b <- glm(survived.ind ~  passengerClass + age + sex -1, 
                  family = binomial(), data = TitanicSurvival)
summary(modmeans_b)

#' 1. Interpret coefficients....
#' 2. Note the number of parameters (5 total)


TitanicSurvival <- TitanicSurvival %>%
  mutate(sex.ind = ifelse(sex == "male", 1, 2),
         class.ind = as.numeric(substr(passengerClass,1,1)))

lrTitanic1<-function(){
  
  # Priors
  
  #for age
  beta ~ dnorm(0, 0.001)
  
  #for sex
  for (i in 1:2) {
    alpha[i] ~ dnorm(0, 0.001)
  }
  
  #for class
  for (i in 1:3) {
    gamma[i] ~ dnorm(0, 0.001)
  }

  
  # Likelihood
  for(i in 1:n){
    logit(p[i]) <-  alpha[sex.ind[i]] + gamma[class.ind[i]] + beta*age.sc[i]
    observed[i]  ~ dbern(p[i])
  }
}


# Bundle data
jagsdata <- list(age.sc = TitanicSurvival$age.sc, 
                 sex.ind = TitanicSurvival$sex.ind, 
                 n = nrow(TitanicSurvival),
                 class.ind = TitanicSurvival$class.ind,
                 observed = TitanicSurvival$survived.ind)


# Parameters to estimate
params <- c("alpha", "beta", "gamma")

# MCMC settings
nc <- 3
ni <- 3000
nb <- 1000
nt <- 2

out.p1 <- jags.parallel(data = jagsdata, parameters.to.save = params, 
                       model.file = lrTitanic1, n.thin= 2, n.chains = 3, 
                       n.burnin = 1000, n.iter = 3000)

MCMCsummary(out.p1)

#' R hat and trace/density plots look ok, but we have an extra parameter!
traplot(out.p1)
denplot(out.p1)
parcorplot(out.p1)  # Oh no!

#' Let's plot the multivariate posterior distribution of say alpha 1 and gamma1
posteriorsamps <- MCMCchains(out.p1)
head(posteriorsamps)
posteriorsamps <- data.frame(posteriorsamps)
ggplot(posteriorsamps, aes(alpha.1., gamma.1.)) + geom_point()
pairs(posteriorsamps)
with(posteriorsamps, plot(alpha.1.[1:50], type = "l"))
with(posteriorsamps, lines(gamma.1.[1:50], type = "l", col = "red"))
     
 
#' ## Session Information
sessionInfo()

 