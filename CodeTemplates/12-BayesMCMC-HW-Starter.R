#' ---
#' title: "Exercise 11.3 Bayesian Inference"
#' output: html_document
#' ---
#' 
#' Load libraries
## ----preamble, warning=FALSE, message=FALSE--------------------------------------
library(R2jags) # for interfacing R and JAGS
library(MCMCvis) # for summarizing JAGS output
library(mcmcplots) # for visualizing MCMC output
library(ggplot2) # for plotting
knitr::opts_chunk$set(fig.width = 6,
               fig.height = 6,
               fig.align = "center")
set.seed(89990)

#' 
#' # Data generation
#' 
#' Let's simulate data using the code from Kery 2010. 
#' 
## ----dataGen---------------------------------------------------------------------
n.groups <- 3 # number of groups
n.sample <- 10 # number of individuals per group
n <- n.groups * n.sample  # Total number of data points
x <- rep(1:n.groups, rep(n.sample, n.groups)) # Indicator for population
pop <- factor(x, labels = c("Pyrenees", "Massif Central", "Jura"))
length <- runif(n, 45, 70)      # Obs. body length (cm) is rarely less than 45

Xmat <- model.matrix(~ pop*length)
beta.vec <- c(-250, 150, 200, 6, -3, -4)

lin.pred <- Xmat[,] %*% beta.vec    # Value of lin.predictor
eps <- rnorm(n = n, mean = 0, sd = 10)  # residuals 
mass <- lin.pred + eps          # response = lin.pred + residual
vip.dat<-data.frame(length = c(length[1:10], length[11:20], length[21:30]), 
                    mass = c(mass[1:10], mass[11:20],mass[21:30]), 
                    pop = pop )
vip.dat

#' 
#' 2. Plot the data:
#' 
## ----fig.height = 6, fig.width = 8, out.width = "65%"----------------------------
ggplot(vip.dat, aes(x = length, y = mass, colour = pop)) + 
  geom_point() + geom_smooth(method = "lm")

#' 
#' 3. The parameters used to simulate the data,
#'  `beta.vec` = (-250, 150, 200, 6, -3, -4), correspond to parameters 
#'   in a model formulated using *effects* coding.  Fit a linear model in
#'   R using *effects* coding (the default) and show that you can recover 
#'   the parameters used to generate the data.
#' 
## --------------------------------------------------------------------------------
lmeffects <- lm(mass ~ pop*length, data = vip.dat)
summary(lmeffects)
ci<-confint(lmeffects)
params <- data.frame(true.betas = beta.vec, CI = ci)
params

#' 
#' All of the 95% confidence intervals include the
#'  values of the parameters used to generate the data.
#' 
#' 4. Fit the same model in R, but using means coding.
#' 
## --------------------------------------------------------------------------------
lmmeans <- lm(mass ~ pop*length -1 -length, data = vip.dat)
summary(lmmeans)

#' 
#' 5. Fit the same model using means coding in JAGS.
#' 
## ----initialJAGS, results="hide"-------------------------------------------------
# Build model
  model <- function(){
    # Priors
      for(i in 1:n.groups){
        alpha[i] ~ dnorm(0, 0.001) # intercepts for each pop
        beta[i] ~ dnorm(0, 0.001) # slopes for each pop
      }
    ## alpha[1] ~ dnorm(0, 0.001)
    ## alpha[2] ~ dnorm(0, 0.001)
    ## alpha[3] ~ dnorm(0, 0.001)
    ## beta[1] ~ dnorm(0, 0.001)
    ## beta[2] ~ dnorm(0, 0.001)
    ## beta[3] ~ dnorm(0, 0.001)
      sigma ~ dunif(0, 100)
      tau <- 1/(sigma*sigma)
    # Likelihood
      for(i in 1:n){
        mass[i] ~ dnorm(mu[i], tau)
        mu[i] <- alpha[popind[i]] + beta[popind[i]]*length[i]
      }
    # Derived quantities
      a.effe2 <- alpha[2] - alpha[1]
      a.effe3 <- alpha[3] - alpha[1]
      b.effe2 <- beta[2] - beta[1]
      b.effe3 <- beta[3] - beta[1]
  }

#' 
#' We need to create a population indicator that goes from 1 to 3, so for each
#' observation, JAGS knows which intercept and slope to use.
## --------------------------------------------------------------------------------
# Bundle Data
  jags.data <- list(mass=as.numeric(mass), popind = as.numeric(pop), 
                    length=length, n.groups=max(as.numeric(pop)), n=n)
  jags.data

#' 
#' Now, fit the model....
#' 
## --------------------------------------------------------------------------------
# Inits Function
  inits <- function(){
    list(alpha=rnorm(n.groups, 0, 2),
         beta=rnorm(n.groups, 1, 1),
         sigma=rlnorm(1))
  }
# Parameters to estimate
  parameters <- c("alpha","beta","sigma",
                  "a.effe2", "a.effe3", "b.effe2", "b.effe3")
# JAGS
  out <- jags(model = model,
              data = jags.data,
              inits = inits,
              parameters.to.save = parameters,
              n.thin = 2,
              n.chains = 3,
              n.iter = 1200,
              n.burnin = 200,
              progress.bar = "none")

#' 
#' 
#' 
## ----initialResults--------------------------------------------------------------
MCMCsummary(out) 

#' 
#' 
#' And now, let's just compare the parameters (from the effects model) to those
#' used to generate the data:
#' 
## --------------------------------------------------------------------------------
beta.vec
MCMCsummary(out, 
            params = c('alpha[1]',  'a.effe2', 'a.effe3', 'beta[1]', 'b.effe2', 'b.effe3'),
            ISB = FALSE, exact = TRUE) 

#' 
