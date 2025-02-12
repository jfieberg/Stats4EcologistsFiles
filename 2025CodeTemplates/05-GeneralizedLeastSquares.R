#' Load libraries
#+ message=FALSE, warning=FALSE
library(knitr)
library(nlme) # for the gls function used to fit gls models
library(dplyr) # for data wrangling
library(ggplot2) # for plotting
library(gridExtra) # for multi-panel plots


#' Jaw data
males<-c(120, 107, 110, 116, 114, 111, 113, 117, 114, 112)
females<-c(110, 111, 107, 108, 110, 105, 107, 106, 111, 111)
jawdat <- data.frame(jaws = c(males, females),
                     sex = c(rep("M",10), rep("F", 10)))


#' Fit GLS model in book 
gls_ttest <- gls(jaws ~ sex, weights = varIdent(form = ~ 1 | sex), data = jawdat)
summary(gls_ttest)


#' Compare to linear model assuming constant variance 
lm_ttest <- gls(jaws ~ sex, data = jawdat)
AIC(lm_ttest, gls_ttest)


#' ## Understanding the variance parameters:
#'
#' log(sigma_s) = log(sigma) + log(delta)*I(sex = female)
#' 
#' log(sigma_s) = gamma0 + gamma1*I(sex=female)
#' 
#' - varpars[1] = log(delta)
#' - varpars[2] = log(sigma)
(varpars<-attr(summary(gls_ttest)$apVar, "Pars"))

#' Sigma for males, females
exp(varpars[2]) # sigma_m
exp(varpars[1]+varpars[2]) #sigma_f


#' Check to see if these make sense 
jawdat %>% group_by(sex) %>%
  dplyr::summarise(sd(jaws))
 

#' ## Sockeye example
library(Data4Ecologists)
data(sockeye)


ggplot(sockeye, aes(MisEsc, SpnEsc, col = Year, shape = Run)) + geom_point() + 
  geom_abline(intercept = 0, slope = 1) +
  labs(x = "Estimate of fish passing Mission", y = "Estimate of spawning escapement + in-river catch")


lmsockeye <- lm(SpnEsc ~ MisEsc, data = sockeye)
ggResidpanel::resid_panel(lmsockeye, plots = c("resid", "qq", "hist"), nrow = 1)

#' Models considered in the book
fixedvar <- gls(SpnEsc ~ MisEsc, weights = varFixed(~ MisEsc), data = sockeye)
varpow <- gls(SpnEsc ~ MisEsc, weights = varPower(form = ~ MisEsc), data = sockeye)
varexp <- gls(SpnEsc ~ MisEsc, weights = varExp(form = ~ MisEsc), data = sockeye) 
varconstp <- gls(SpnEsc ~ MisEsc, weights = varConstPower(form = ~ MisEsc), data = sockeye)

#' Also fit a standard linear model and compare
lmsockeye <- gls(SpnEsc ~ MisEsc, data = sockeye)
AIC(lmsockeye, fixedvar, varpow, varexp, varconstp)

#' Lets look at the model with smallest AIC
summary(varconstp)

#' Variance for a particular observation
sig2i <- 0.172081*(119.49037 + abs(sockeye$MisEsc)^1.10369)

#' Standardized residuals should have constant variance
sockeye <- sockeye %>%
  mutate(fitted = predict(varconstp),
         stdresids  = (SpnEsc - predict(varconstp))/sig2i)

#' Let's see if that is the case
ggplot(sockeye, aes(fitted, stdresids)) + geom_point() +
  geom_hline(yintercept = 0) 


#' There is also a default plot that does this for us
plot(varconstp)

#' We can also see if they look Normally distributed
ggplot(sockeye, aes(sample = stdresids)) +
  stat_qq() +
  stat_qq_line()

#' Scale-location plot
sockeye <- sockeye %>% 
  mutate(absstdresids  = sqrt(abs(stdresids)))

ggplot(sockeye, aes(fitted, absstdresids)) + geom_point() +
  geom_smooth() 

#' More complicated models can also be considered
varconstp2 <- gls(SpnEsc ~ MisEsc, weights = varConstPower(form = ~ MisEsc | Run), 
                  data = sockeye)
summary(varconstp2)
AIC(varconstp, varconstp2)


#' Model where the variance depends on the mean
varmean <- gls(SpnEsc ~ MisEsc, weights = varPower(form = ~ fitted(.)), 
                  data = subset(sockeye, Run=="ESum")) 
summary(varmean)



# Design matrix for our observations
xmat <- model.matrix(~ MisEsc, data=sockeye)

# Regression coefficients
betahat<-coef(varconstp)

# Predictions
sockeye$SpnEsc.hat <- predict(varconstp)
cbind(head(xmat%*%betahat), head(sockeye$SpnEsc.hat))


# Sigma^
Sigmahat <- vcov(varconstp)

# var/cov(beta0 + beta1*X)
varcovEYhat<-xmat%*%Sigmahat%*%t(xmat)

# Pull off the diagonal elements and take their sqrt to 
# get SEs that quantify uncertainty associated with the line
SEline <- sqrt(diag(varcovEYhat))


# Confidence interval for the mean
sockeye$upconf <- sockeye$SpnEsc.hat + 1.645 * SEline
sockeye$lowconf <- sockeye$SpnEsc.hat - 1.645 * SEline
p1 <- ggplot(sockeye, aes(MisEsc, SpnEsc)) + geom_point() +
  geom_line(aes(MisEsc, SpnEsc.hat)) +
  geom_line(aes(MisEsc, lowconf), col = "red", lty = 2) +
  geom_line(aes(MisEsc, upconf), col = "red", lty = 2) +
  xlab("Estimate at Mission") + ylab("Estimate of C+Spawners")
p1


# Approx. prediction interval, ignoring the uncertainty in theta
varpredict <- SEline ^ 2 + sig2i ^ 2
sockeye$pupconf <- sockeye$SpnEsc.hat + 1.645 * sqrt(varpredict)
sockeye$plowconf <- sockeye$SpnEsc.hat - 1.645 * sqrt(varpredict)
ggplot(sockeye, aes(MisEsc, SpnEsc)) + geom_point() +
  geom_line(aes(MisEsc, SpnEsc.hat)) +
  geom_line(aes(MisEsc, lowconf), col = "red", lty = 2) +
  geom_line(aes(MisEsc, upconf), col = "red", lty = 2) +
  geom_line(aes(MisEsc, plowconf), col = "blue", lty = 2) +
  geom_line(aes(MisEsc, pupconf), col = "blue", lty = 2) +
  xlab("Estimate at Mission") + ylab("Estimate of C+Spawners")

