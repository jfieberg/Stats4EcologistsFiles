#' # Bayesian Linear Model

#' ## Preamble
#' 
#' Load packages
#+warning=FALSE, message=FALSE
library(R2jags) # for interfacing with R
library(MCMCvis) # for summarizing MCMC output
library(mcmcplots) # for plotting MCMC output
library(ggplot2) # for other plots
library(patchwork) # for multi-panel plots
library(dplyr) # for data wrangling
library(WVPlots) # for shading of posterior density

#' Write out the model
linreg<-function(){

# Priors
 beta0 ~ dnorm(0,0.001)
 beta1 ~ dnorm(0,0.001)
 sigma ~ dunif(0, 100)

# Derived quantities
 tau <- 1/ (sigma * sigma)

# Likelihood
 for (i in 1:n) {
    age[i] ~ dnorm(mu[i], tau) 
    mu[i] <- beta0 + beta1*proportion.black[i]
 }

# Assess model fit using a sums-of-squares-type discrepancy
 for (i in 1:n) {
    residual[i] <- age[i]-mu[i]       # Residuals for observed data
    sq[i] <- pow(residual[i], 2)    # Squared residuals for observed data

# Generate replicate data and compute fit stats for them
    y.new[i] ~ dnorm(mu[i], tau) # one new data set at each MCMC iteration
    sq.new[i] <- pow(y.new[i]-mu[i], 2)  # Squared residuals for new data
 }
 fit <- sum(sq[])           # Sum of squared residuals for actual data set
 fit.new <- sum(sq.new[])       # Sum of squared residuals for new data set
}

#' Read in data and gather everything needed to fit the model in JAGS
#+warning=FALSE, message = FALSE
library(abd)
data("LionNoses")

# Bundle data
jags.data <- list(age = LionNoses$age, 
                  proportion.black = LionNoses$proportion.black, 
                  n = nrow(LionNoses))

#' Parameters to estimate
params <- c("beta0", "beta1", "sigma", "mu", "y.new", "fit", "fit.new", "residual")

# MCMC settings
nc = 3
ni = 2200
nb = 200
nt = 1

# Start Gibbs sampler
lmbayes <- jags(data = jags.data,  
                parameters.to.save = params, 
                model.file = linreg, 
                n.thin = nt, 
                n.chains = nc, 
                n.burnin = nb, 
                n.iter = ni, 
                progress.bar="none")

#' Make sure we have reached convergence and density and traceplots look OK
MCMCsummary(lmbayes, params = c("beta0", "beta1", "sigma"))

#+ fig.height=8, fig.width=4
denplot(lmbayes, parms = c("beta0", "beta1", "sigma"), ci=0.95)
traplot(lmbayes, parms = c("beta0", "beta1", "sigma"))

#' Fit the same model using lm and compare results
lmnoses <- lm(age ~ proportion.black, data = LionNoses)
summary(lmnoses)
MCMCsummary(lmbayes, params = c("beta0", "beta1", "sigma"))


#' Residual plots using posterior residuals and fitted values
lmresids <- MCMCpstr(lmbayes, params = "residual", func = mean)
lmfitted <- MCMCpstr(lmbayes, params = "mu", func = mean)
jagslmfit <- data.frame(resid = lmresids$residual, fitted = lmfitted$mu)
jagslmfit$std.abs.resid <- sqrt(abs(jagslmfit$resid/sd(jagslmfit$resid)))


#' Compare residual plots 
#+bayesresids, fig.height=4, fig.width=12, out.width = "100%", fig.align='center'
p1 <- ggplot(jagslmfit, aes(fitted, resid)) + geom_point() + 
         geom_hline(yintercept = 0) + geom_smooth()
p2 <- ggplot(jagslmfit, aes(sample = resid)) + stat_qq() + stat_qq_line()
p3 <- ggplot(jagslmfit, aes(fitted, std.abs.resid)) + geom_point() + 
         geom_smooth() + ylab("sqrt(|Standardized Residuals|)")
p1 + p2 + p3


#+freqresids, fig.height=4, fig.width=12, out.width = "100%"
ggResidpanel::resid_panel(lmnoses, plots = c("resid", "qq", "ls"), nrow = 1, smoother = TRUE)


#' Goodness of fit test 
fitstats <- MCMCpstr(lmbayes, params = c("fit", "fit.new"), type = "chains") 
T.extreme <- fitstats$fit.new >= fitstats$fit
(p.val <- mean(T.extreme))

#' Posterior for the difference in fits
#+out.width="60%", fig.width=8, fig.height=5
fitstats <- MCMCchains(lmbayes, params = c("fit", "fit.new")) 
fitstats <- as.data.frame(fitstats) %>% mutate(postdiff = fit.new-fit)
WVPlots::ShadedDensity(frame = fitstats, 
                       xvar = "postdiff",
                       threshold = 0,
                       title = "Posterior distribution: SSE(sim data)-SSE(obs data)",
                       tail = "right")+
  annotate("text", x=85, y = 0.005, 
           label="Better fit to \n observed data")+
  annotate("text", x=-40, y = 0.005, 
           label="Better fit to \n simulated data") 

#' ## Confidence and prediction intervals
#'
#' ### Confidence intervals
#' Store the MCMC samples from the posterior.  This can be done with
#' MCMCpstr (which stores the result as a list) or MCMCchains (which 
#' stores the results as a matrix). 
betas <- MCMCpstr(lmbayes, params = c("beta0", "beta1"), type = "chains")

#' Create a vector of observations where we want to calculate CIs and PIs 
prop.black <-seq(from = min(LionNoses$proportion.black), 
                 to = max(LionNoses$proportion.black), 
                 length = 10)

# Number of MCMC samples
nmcmc <- dim(betas$beta0)[2]

# Number of unique values of proportion.black where we want predictions
nlions <- length(prop.black)

# Set up a matrix to hold 95% credible intervals for each value of prop.black
conf.int <- matrix(NA, nlions, 2) 

# Loop over values for proportion.black
for(i in 1:length(prop.black)){
  # Estimate the mean age for the current value of prop.black and for 
  #   each MCMC sample of beta0 and beta1
  age.hats <- betas$beta0 + rep(prop.black[i], nmcmc)*betas$beta1
  conf.int[i,] <- quantile(age.hats, prob = c(0.025, 0.975))
}


#+ fig.align='center', out.width = "55%", fig.height = 5, fig.width =6
beta.hat <- MCMCpstr(lmbayes, params = c("beta0", "beta1"), func = mean)
phats <- data.frame(est = rep(beta.hat$beta0, nlions) + rep(beta.hat$beta1, nlions)*prop.black, 
                    LCL = conf.int[,1], 
                    UCL = conf.int[,2], proportion.black = prop.black)
ggplot(phats, aes(proportion.black, est)) +
  geom_ribbon(aes(ymin = LCL, ymax = UCL), fill = "grey70", alpha = 2) +
  geom_line() +ylab("Age") + 
  geom_point(data = LionNoses, aes(proportion.black, age))


#' ### Prediction intervals
#' 
# Pull off the samples from the posterior distribution of sigma
sigmap <- MCMCpstr(lmbayes, params = "sigma", type = "chains")$sigma

# Set up a matrix to hold 95% prediction intervals for each value of prop.black
pred.int <- matrix(NA, nlions, 2) 

# Loop over values for proportion.black
for(i in 1:length(prop.black)){
  # Estimate the mean age for the current value of prop.black and for
  #    each MCMC sample of beta0 and beta1
  age.inds <- betas$beta0 + rep(prop.black[i], nmcmc)*betas$beta1 + 
                     rnorm(nmcmc, mean = 0, sd = sigmap)
  # Now for prediction intervals using the quantiles of the age values
  pred.int[i,] <- quantile(age.inds, prob = c(0.025, 0.975))
}


#+ bayesplotcred1, fig.align='center', out.width = "55%", fig.height = 5, fig.width =6 
beta.hat <- MCMCpstr(lmbayes, params = c("beta0", "beta1"), func = mean)
phats <- data.frame(phats, 
                    LPL = pred.int[,1], 
                    UPL = pred.int[,2])
ggplot(phats, aes(proportion.black, est)) +
  geom_line(aes(proportion.black, LCL), col = "red", lty = 2) +
  geom_line(aes(proportion.black, UCL), col = "red", lty = 2) +
  geom_line(aes(proportion.black, LPL), col = "blue", lty = 2) +
  geom_line(aes(proportion.black, UPL), col = "blue", lty = 2) +
  geom_line() +ylab("Age") + 
  geom_point(data = LionNoses, aes(proportion.black, age))



#' ## CIs and PI's the easier way 
CI <- MCMCpstr(lmbayes, params = c("mu"), 
               func = function(x) quantile(x, probs = c(0.025, 0.975)))$mu
PI <- MCMCpstr(lmbayes, params = c("y.new"), 
               func = function(x) quantile(x, probs = c(0.025, 0.975)))$y.new



#+bayesplotcred2, fig.align='center', out.width = "55%", fig.height = 5, fig.width =6 
phats2 <- data.frame(est =MCMCpstr(lmbayes, params = c("mu"), func = mean)$mu, 
                     proportion.black = LionNoses$proportion.black,
                    LCL = CI[,1], 
                    UCL = CI[,2], 
                    LPL = PI[,1], 
                    UPL = PI[,2])
ggplot(phats2, aes(proportion.black, est)) +
  geom_line(aes(proportion.black, LCL), col = "red", lty = 2) +
  geom_line(aes(proportion.black, UCL), col = "red", lty = 2) +
  geom_line(aes(proportion.black, LPL), col = "blue", lty = 2) +
  geom_line(aes(proportion.black, UPL), col = "blue", lty = 2) +
  geom_line() +ylab("Age") + 
  geom_point(data = LionNoses, aes(proportion.black, age))

