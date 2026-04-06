#' ---
#' title: "Logistic regression models"
#' author: "John Fieberg"
#' output: 
#'    html_document:
#'      toc: true
#'      toc_depth: 2
#'      toc_float: true
#' ---
#' 
## ----echo=FALSE, warning=FALSE,message=FALSE--------------------
set.seed(12111)
library(kableExtra) # for tables
library(patchwork) # for multi-panel plots
library(ggplot2) # for plotting
library(dplyr) # for data wrangling 
library(car) # for residual plots
library(modelsummary) # for summarizing models
library(MASS) # for generating multivariate random normal variables
library(ggeffects) # for summarizing fitted models 
library(effects) # for summarizing fitted models 
library(performance) # for residual plots
library(ggthemes) # for colorblind pallete
knitr::opts_chunk$set(fig.align = 'center',
                      out.width = "55%",
                      fig.height = 5,
                      fig.width = 6)

#' Read in a moose data set  
library(SightabilityModel)
data(exp.m)
str(exp.m)


#' ## Model 1: logit(p) = beta0 + beta1*voc
#'
mod1 <- glm(observed ~ voc, data = exp.m, family = binomial())
summary(mod1)

#' ### Predictions and confidence intervals
#' 
#' Let's begin by considering our initial model containing only `voc`.  
#' We can summarize models in terms of:
#' 
#' $\hat{E}[Y_i|X_i] =\hat{p}_i$ for a range of `voc` values. 
#' 
# Predictions link scale
voc.p<-seq(1, 99, length = 100)
newdat = data.frame(voc = voc.p)
phat1<-predict(mod1, type = "link", se.fit = T, newdata = newdat)
lcl.l <- plogis(phat1$fit + 1.96 * phat1$se.fit)
ucl.l <- plogis(phat1$fit - 1.96 * phat1$se.fit)
pe.l <- plogis(phat1$fit)

#' Predictions on response scale using delta method for SE
phat2<-predict(mod1, type = "response", se.fit = T, newdata = newdat)
lcl.r <- phat2$fit + 1.96 * phat2$se.fit
ucl.r <- phat2$fit - 1.96 * phat2$se.fit
pe.r <- phat2$fit

# Combine and plot
#+ fig.alt = "Plot of prob(detect) with confidence intervals
newdat <- cbind(newdat, phat = pe.r, lcl.r, ucl.r, phat2 = phat2$fit, lcl.l, ucl.l, pe.l)
ggplot(newdat, aes(voc, phat)) + geom_line() + 
  geom_line(aes(voc, lcl.l, color = "logit scale"), lty = 2) +
  geom_line(aes(voc, ucl.l, color = "logit scale"), lty = 2) +
  geom_line(aes(voc, lcl.r,color = "response scale"), lty = 2) +
  geom_line(aes(voc, ucl.r, color = "response scale"), lty = 2) + 
  scale_color_manual(name = "CI type", 
                     values = c("logit scale" = "blue", "response scale" = "red")) +
  xlab("Visual Obstruction") + ylab(expression(hat(p)))


#' ## Effect plots: Visualizing generalized linear models 
#'  
#' Let's begin by creating our own plots, which will give us insights 
#' into how the effects/ggeffects packages work using a slightly more complicated 
#' model.

#' Fit model with year (as factor) and voc
exp.m$year <- as.factor(exp.m$year)
mod2 <- glm(observed ~ voc + year, data = exp.m, family = binomial())
summary(mod2)


#'
#' ## Logit versus Probability scale 
#'  
#' Note:  See Figure 16:9 in the book. Effects are additive on the logit 
#' scale but not on the probability scale.
newdat<-expand.grid(voc=seq(0,100,5), year=unique(exp.m$year))  
newdat$p.hats<-predict(mod2, newdat=newdat, type="response")
newdat$logit.p.hats<-predict(mod2, newdata=newdat, type="link")
plot.logitp<-ggplot(newdat, aes(voc, logit.p.hats, colour=year))+geom_line()+
  ylab("logit(p)")
plot.p<-ggplot(newdat, aes(voc, p.hats, colour=year))+geom_line()+ylab("p")
plot.logitp + plot.p 

#' 
#' 
#' ### Effect plots 
#' 
#' Vary 1 predictor while holding all other predictors at a common value
#' (e.g., the mean for continuous variables and the modal value for
#' categorical variables). 
#' 
#' The `effect` function in the effects` package:
#' 
#' - Fixes all continuous variables (other than the one of interest)
#'  at their mean values
#' - For categorical predictors, it averages predictions on the link scale,
#'  weighted by the proportion of observations in each category
#' 
#' The effect plot for year depicts:  
#'
#'   P(detect | year[varied], voc = mean(voc)).
#' 
#' The effect plot for voc depicts:  
#' 
#'   sum(over years, P(detect | voc[varied], year = year[i])*(n year i)/(total n)
#' 
ggeffect(mod2, "year") # Prediction for each year, with voc set to its mean value
ggeffect(mod2, "voc") # Prediction for a range of vocs, averaged over years

#'  
#+ out.width = "100%", fig.height = 4, fig.width = 8 
p1 <- plot(ggeffect(mod2, "year"))
p2 <- plot(ggeffect(mod2, "voc"))
p1 + p2



#' 
#' ###  Predictions using the `ggpredict` function
#' 
#' Whereas the `effects` function averages predictions across
#'  the different levels of a categorical predictor, producing what 
#'  are sometimes referred to as *marginal effects*, the `ggpredict` 
#'  function will provide *adjusted predictions* where all variables
#'  except a focal variable remain fixed at their mean (or median), modal, or 
#'  user-specified values.  
ggpredict(mod2, "voc") # Prediction for range of voc's with year set to 2005
ggpredict(mod2, "year") # Prediction for each year with voc set to its MEDIAN value
median(exp.m$voc)
mean(exp.m$voc)

#' To match up ggpredict and ggeffect for year:
ggpredict(mod2, terms = "year", condition = c(voc=mean(exp.m$voc)))
ggeffect(mod2, "year")

#' Or, we can predict for a range of voc's for each year by specifying
#' both terms.
#+ out.width = "100%", fig.height = 4, fig.width = 12 
adjustpredict <- ggpredict(mod2, terms = c("voc", "year"))
plot(adjustpredict, grid = TRUE)  
ggpredict(mod2, terms = c("voc", "year"))

