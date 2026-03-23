#' HW 4
library(Data4Ecologists)
data(beepollen)

#' Load packages
#+ warning = FALSE, message = FALSE
library(emdbook) # for delta method function
library(dplyr) # for data wrangling
library(ggplot2) # for plotting


#' Lets fit 2 models, one in which the only the asympotote differs and one
#' where the rate of increase also differs.
neg.likebee1 <- function(params, data){
  beta1 <- params[1]
  beta2 <- params[2]
  beta3 <- params[3]
  logsigma <- params[4]
  sigma<-exp(logsigma)
  EY <- (beta1 + beta2*data$queen)*(1- exp(-beta3*data$duration))
  -sum(dnorm(data$removed, mean = EY, sd = sigma, log = TRUE))
}


neg.likebee2 <- function(params, data){
  beta1 <- params[1]
  beta2 <- params[2]
  beta3 <- params[3]
  beta4 <- params[4]
  logsigma <- params[5]
  sigma<-exp(logsigma)
  EY <- (beta1 + beta2*data$queen)*(1- exp(-(beta3+beta4*data$queen)*data$duration))
  -sum(dnorm(data$removed, mean = EY, sd = sigma, log = TRUE))
}

fitbee1 <- optim(par=c(0.75, 0, 0.1, log(0.2)), 
               fn = neg.likebee1, dat = beepollen, 
               method="BFGS", hessian=TRUE)

fitbee2 <- optim(par=c(0.75, 0, 0.1, 0, log(0.2)), 
                 fn = neg.likebee2, dat = beepollen, 
                 method="BFGS", hessian=TRUE)

#' Calculate AIC's = -2*Log Likelihood + 2*number of parameters
2*fitbee1$value + 2*length(fitbee1$par)
2*fitbee2$value + 2*length(fitbee2$par)

#' Plot the data with both models overlayed
betas1 <- fitbee1$par # optim parameters
betas2 <- fitbee2$par

newdat <- data.frame(expand.grid(duration = seq(min(beepollen$duration), max(beepollen$duration), length = 100), 
                                 queen = c(0, 1)))
newdat$fit1 <- (betas1[1] + betas1[2]*newdat$queen)*(1 - exp(-betas1[3]*newdat$duration))
newdat$fit2 <- (betas2[1] + betas2[2]*newdat$queen)*(1 - exp(-(betas2[3]+betas2[4]*newdat$queen)*newdat$duration))

ggplot(beepollen, aes(duration, removed)) + geom_point() + 
  geom_line(data = newdat, aes(duration, fit1), col="red") + 
  geom_line(data = newdat, aes(duration, fit2), col="blue") +
  facet_wrap(~queen)
