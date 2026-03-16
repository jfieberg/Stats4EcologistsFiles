#' # JAGS t-test/linear model

#' ## Preamble

#' To install BEST package (since no longer on CRAN)
#remotes::install_version("BEST", version = "0.1.0")

#' Increase width of output
options(width=150)

#' Install & Load any libraries  
#+ warning=FALSE, message=FALSE
library(rjags) # for accessing JAGS
library(R2jags) # for interfacing with JAGS
#library(mcmcplots) # for plotting results (no longer on CRAN)
library(ggplot2) # for some later plotting
library(MCMCvis) # for summarizing MCMC results, other options: ggmcmc, bayesplot and basicMCMCplots
#library(BEST) # for prior / posterior overlap (no longer on CRAN)
library(dplyr) # for data wrangling
library(knitr) # for options when knitting
library(ggplot2) # for plotting
library(ggthemes) # for colorblind theme
 
#' Read in data
males<-c(120, 107, 110, 116, 114, 111, 113, 117, 114, 112)
females<-c(110, 111, 107, 108, 110, 105, 107, 106, 111, 111)


#' JAGS/Bugs model for jaw data
jaw.mod<-function(){
  
  # Likelihood (Y | mu) = normal(mu[sex], sigma^2)  
    for(i in 1:nmales){
      males[i] ~ dnorm(mu.male, tau) 
    }
    for(i in 1:nfemales){
      females[i] ~ dnorm(mu.female, tau)
    }
  
  # Priors 
  mu.male ~ dnorm(100, 0.001) # mean of male jaw lengths
  mu.female ~ dnorm(100, 0.001) # mean of female jaw lengths
  sigma ~ dunif(0, 30) # common sigma
  tau <- 1/(sigma*sigma) #precision
  
  # Derived quantities:  difference in means
    mu.diff <- mu.male - mu.female
}
  

# Function to generate initial values
init.vals<-function(){
    mu.male <- rnorm(1, 100, 100)
    mu.female <- rnorm(1, 100, 100)
    sigma <- runif(1, 0, 10) 
    out <- list(mu.male = mu.male, 
                mu.female = mu.female, 
                sigma = sigma)
}


#' Create rest of the data for the model
nmales<-length(males)
nfemales<-length(females)

 

#' If spinning, add "progress.bar="gui" or "none" to jags function
#' (for Mac users, use the "none" option!)
#+ message=FALSE, warning=FALSE, results="hide"   
t.test.jags <- jags(data=c("males", "females","nmales",  "nfemales"), 
                    inits=init.vals,
                    parameters.to.save=c("mu.male", "mu.female","sigma", "mu.diff"),
                    progress.bar = "none",
                    n.iter=10000, 
                    n.burnin=5000,
                    n.thin=1, 
                    n.chains=3,
                    model.file=jaw.mod)

#' Or, to take advantage of parallel processing
t.test.jags<-jags.parallel(data=c("males", "females","nmales",  "nfemales"), 
                           inits=init.vals,
                           parameters.to.save=c("mu.male", "mu.female","sigma", "mu.diff"),
                           n.iter=10000, 
                           n.burnin=5000,
                           n.thin=1, 
                           n.chains=3,
                           model.file=jaw.mod)

#' Look at the output
t.test.jags


#" Summarize output for a subset of parameters
MCMCsummary(t.test.jags, params = c("mu.male", "mu.female"))


# Use MCMCsummary to pull off posterior means
bayesests <- MCMCpstr(t.test.jags, params = c("mu.male", "mu.female"), func = mean)
bayesests

# Use MCMCsummary to pull of upper and lower 95% credible interval limits
bayesci <-  MCMCpstr(t.test.jags, params = c("mu.male", "mu.female"), 
                        func = function(x) quantile(x, probs = c(0.025, 0.975)))
bayesci

#' Look at posterior density and trace plots. ind = TRUE plots all 3 chains.
#+ fig.alt = "Density and trace plots"
MCMCtrace(t.test.jags, ind = TRUE, , pdf = FALSE)

#' Caterplot to show 68% and 95% credible intervals
MCMCplot(t.test.jags, 
          params=c("mu.male", "mu.female")) 

#' # Compare to linear model fit using lm

# Fit linear model in frequentist framework
jawdat<-data.frame(jaws=c(males, females), 
                    sex=c(rep("Male", length(males)), 
                          rep("Female", length(females))))
lm.jaw<-lm(jaws~sex-1, data=jawdat)


# Store results
betaf <- coef(lm.jaw)
CIf <-confint(lm.jaw) 
ests<-data.frame(estimate = c(bayesests$mu.female, bayesests$mu.male, betaf), 
                 LCL = c(bayesci$mu.female[1],   bayesci$mu.male[1], CIf[,1]), 
                 UCL = c(bayesci$mu.female[2],   bayesci$mu.male[2], CIf[,2]), 
                 param = c("Mean females", "Mean males"),
                 Method = rep(c("Bayesian", "Frequentist"), each = 2))


#' Compare Bayesian and frequentist estimates
ggplot(ests, aes(param,estimate, col = Method)) + 
  geom_point(position = position_dodge(width = 0.2))+ 
  geom_pointrange(aes(ymin = LCL, ymax= UCL), position = position_dodge(width = 0.2))+
  ylab("Estimate") + xlab("") +
  scale_x_discrete(labels = c('Mean females' = expression(mu[f]),
                              'Mean males'   = expression(mu[m]))) + 
  scale_colour_colorblind()+
  theme(text = element_text(size = 20))

#' For more control, we can also access the samples for the different chains using MCMCchains
mus<- MCMCchains(t.test.jags, params=c("mu.male", "mu.female")) 
str(mus)
head(mus)

#' Plot posteriors
#+ fig.alt = "Posterior densities for mean jaw length"
ggplot(mus) +
  geom_density(aes(mu.male, col = "Male")) +
  geom_density(aes(mu.female, col = "Female")) +
  scale_color_manual(values = c("Male" = "red", "Female" = "blue")) +
  labs(color = "Sex", x = expression(mu), y = "Density")
 

#' Another option
#+ fig.alt = "Posterior density plots"
musMCMC2<-data.frame(mu = c(mus[,1], mus[,2]),
                            sex=rep(c("Male", "Female"), each = 15000))
ggplot(musMCMC2, aes(mu, col=sex))+geom_density() 


#' ## Prior posterior overlap
#' 
#' We can look at prior/posterior overlap using `MCMCtrace`. This requires that
#' we first generate samples from our priors, ordered alphabetically in a matrix.
#' We also need the same number of iterations as captured by the posterior distribution.


#' Generate the same number of samples from the prior distribution for sigma
niters<-15000 # total number of iterations
Prior.samples<-runif(niters, 0, 30)  

#' Now, call MCMCtrace
MCMCtrace(t.test.jags, params = c("sigma"),
          priors = Prior.samples,
          pdf = FALSE,
          post_zm = FALSE)

#' Or, create our own plot by combining prior and posterior samples in a data frame
#+ fig.alt = "Prior and posterior distributions for sigma"
priorpost <- data.frame(sigma = c(MCMCchains(t.test.jags, params="sigma"), Prior.samples),
                        distribution = rep(c("Posterior", "Prior"), each = 15000))

ggplot(priorpost, aes(sigma, col=distribution)) +
   geom_density()                  
 
