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
 
#+ fig.align="center",  fig.height = 5, fig.width = 6, out.width = "55%"
ggplot(exp.m, aes(voc,observed))+theme_bw()+
    geom_point(position = position_jitter(w = 2, h = 0.05), size=3) +
    geom_smooth(colour="red") + xlab("Visual Obstruction") +
    ylab("Detection = 1")

#' Fit logistic regression model with voc
mod1 <- glm(observed ~ voc, data = exp.m, family = binomial())
summary(mod1)

#' Inference based assumed Normal distribution for the sampling distribution
beta1 <- coef(mod1)[2]
SEbeta1 <- sqrt(diag(vcov(mod1)))[2]
oddsr <- exp(beta1)
CI_beta <- beta1 + c(-1.96, 1.96)*SEbeta1
CI_odds <- exp(CI_beta)
 
#' Alternatively, we could calculate profile-likelihood based
#'  confidence intervals by inverting the likelihood ratio test 
#'  Recall, profile-likelihood intervals include all values, $\widetilde{\beta}$
#'  for which we would **not** reject the null hypothesis 
#'  $H_0: \beta=\widetilde{\beta}$. In small samples, profile likelihood
#'  intervals should perform better than the Normal-based intervals, 
#'  and thus, this is the default approach used by the function `confint`: 
#' 
(ci.prof <- confint(mod1))
(exp(ci.prof)) 

#' Lastly, we can visualize the model using ``ggplot` using `stat_smooth`
ggplot(exp.m, aes(voc,observed)) + theme_bw() + 
    geom_point(position = position_jitter(w = 2, h = 0.05), size=3) +
    xlab("Visual Obstruction") + 
    stat_smooth(method="glm", method.args = list(family = "binomial") )+
    ylab("Detection = 1") 

#' ### Model with continuous and categorical variables 
#' 
#' Lets now consider a second model, where we also include year of the
#' observation as a predictor. 
exp.m$year <- as.factor(exp.m$year)
mod2 <- glm(observed ~ voc + year, data = exp.m, family = binomial())
summary(mod2)


 
#' ### Interaction model
#' 
#' What if we were to include an interaction between `voc` and `year`?
mod3 <- glm(observed ~ voc * year, data = exp.m, family = binomial())
summary(mod3)

#' Write out the design matrix for these observations
exp.m[c(1,37,53,74), c("voc", "year")]

#' Check your understanding:
model.matrix(mod3)[c(1,37,53,74),]

 
#' ## Evaluating assumptions and fit 
#' 
#' ### Residual plots 
#'
#' Pearson residuals versus the linear predictor 
residualPlot(mod1)

#' Bin values of the linear predictor and residuals and then plot
binplot<-binned_residuals(mod1)
plot(binplot)


#' ### Goodness-of-fit test
#'
#' Can again use the sum of Pearson residuals for a GOF test:
#'  
#' $$\chi^2_{n-p} = \sum_{i=1}^n\frac{(Y_i-E[Y_i|X_i])^2}{Var[Y_i|X_i]}$$
#' 
#' For binary data analyzed using logistic regression, we have:
#' 
#' - $E[Y_i|X_i] = p_i$ 
#' - $Var[Y_i|X_i] =$  $p_i(1-p_i)$
#' 
#' Set up simulation
nsims<-10000
gfit.sim<-gfit.obs<-matrix(NA, nsims, 1)
nobs<-nrow(exp.m)

#' Generate coefficients from estimated sampling distribution so that
#' we can average over uncertaining
beta.hat<-mvrnorm(nsims, coef(mod1), vcov(mod1))
xmat<-model.matrix(mod1)
for(i in 1:nsims){

# Estimate p for each observation
    ps<- 
    
# Generate new random values using the assumed model    
    new.y<-

#' Calculate sum of Pearson residuals        
    gfit.sim[i,]<-sum((new.y-ps)^2/ 
    gfit.obs[i,]<-sum((exp.m$observed-ps)^2/ 
  }
    mean(gfit.sim > gfit.obs) 

#' 
#' ### Aside: alternative methods for goodness-of-fit testing 
#' 

#' Similar test, but using binned data (Hosmer-Lemeshow test). There 
#' are several packages that implement this test. Results will depend somewhat
#' on the number of bins you choose.
ResourceSelection::hoslem.test(exp.m$observed, fitted(mod1), g = 8)
performance::performance_hosmer(mod1)
performance::performance_hosmer(mod1, n_bins=8)
 
#' ## Model comparisons   
#' 
#' Nested models: use likelihood ratio tests based on $\Chi^2$ distribution
drop1(mod2, test="Chisq")
car::Anova(mod2) 
AIC(mod1, mod2, mod3)

#' The model with `voc` and `year` but not their interaction has the lowest AIC.
#'



#' ## Effect plots: Visualizing generalized linear models 
#'  
#' Let's begin by creating our own plots, which will give us insights 
#' into how the effects/ggeffects packages work.
#' 
#' ### Predictions and confidence intervals
#' 
#' Let's begin by considering our initial model containing only `voc`.  
#' We can summarize models in terms of:
#' 
#' $\hat{E}[Y_i|X_i] =\hat{p}_i$ for a range of `voc` values. 
#' 
voc.p<-seq(1, 99, length = 100)
p.hat <- predict(mod1, newdata = data.frame(voc = voc.p), type = "response")
plot(voc.p, p.hat, xlab = "VOC", ylab = "Pr(detect | voc)", type = "l")

#' Now, with CI's:
#' 
# Predictions link scale
phat2<-predict(mod1, type = "link", se.fit = T, newdata = newdat)
lcl.l <- plogis(phat2$fit + 1.96 * phat2$se.fit)
ucl.l <- plogis(phat2$fit - 1.96 * phat2$se.fit)
pe.l <- plogis(phat2$fit)

# Combine and plot
newdat <- cbind(newdat, phat = phat$fit, lcl.r, ucl.r, phat2 = phat2$fit, lcl.l, ucl.l, pe.l)
ggplot(newdat, aes(voc, phat)) + geom_line() + 
  geom_line(aes(voc, lcl.l), lty = 2, col = "blue") +
  geom_line(aes(voc, ucl.l), lty = 2, col = "blue") +
  geom_line(aes(voc, lcl.r), lty = 2, col = "red") +
  geom_line(aes(voc, ucl.r), lty = 2, col = "red") + 
  xlab("Visual Obstruction") + ylab(expression(hat(p)))

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
#'  except a focal variable remain fixed at their mean, modal, or 
#'  user-specified values.  
ggpredict(mod2, "year")
ggpredict(mod2, "voc")

#' To match up ggpredict and ggeffect for year:
ggpredict(mod2, terms = "year", condition = c(voc=mean(exp.m$voc)))
ggeffect(mod3, "year")

#' Or, we can predict for a range of voc's for each year by specifying
#' both terms.
#+ out.width = "100%", fig.height = 4, fig.width = 12 
adjustpredict <- ggpredict(mod2, terms = c("voc", "year"))
plot(adjustpredict, grid = TRUE)  

 
 