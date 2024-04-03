#' ---
#' title: "17-LinearMixedEffectsModels2.R"
#' author: "John Fieberg"
#' output: 
#'    html_document:
#'      toc: true
#'      toc_depth: 2
#'      toc_float: true
#' ---

 
#' ## Load libraries
#+ warning=FALSE,message=FALSE
library(ggplot2);theme_set(theme_bw())
knitr::opts_chunk$set(fig.align='center', 
                      out.width = "55%", 
                      fig.height = 5, 
                      fig.width =6,
                      error = FALSE)
library(tidyverse) # for data wrangling
library(gridExtra) # for multi-panel plots
library(lme4)  # for fitting random-effects models
library(nlme) # for fitting random-effects models
library(glmmTMB) # for fitting random-effects models
library(sjPlot) # for visualizing fitted models
library(modelsummary) # for creating output tables
library(kableExtra) # for tables 
options(kableExtra.html.bsTable = T)
library(ggplot2)# for plotting
library(performance) # for model diagnostics

#' Load data sets
#+warning=FALSE, message=FALSE 
library(Data4Ecologists) # for data
data("RIKZdat")


#' Note that NAP is a "level-1" covariate (it varies within each cluster),
#' whereas exposure is a level-2 covariate. 
boxplot(NAP~Beach, data=RIKZdat, xlab="Beach", ylab="NAP")	
xtabs(~ exposure + Beach, data=RIKZdat)	

#' Only 1 beach has the lowest exposure level:  Recode exposure to have only 2 levels.	
RIKZdat$exposure.c<-"High"	  
RIKZdat$exposure.c[RIKZdat$exposure%in%c(8,10)]<-"Low"	

#' Plot data again
#+ out.width="90%"
ggplot(RIKZdat, aes(NAP, Richness, col=exposure.c)) +
  geom_point()+
  geom_smooth(method="lm", formula = y ~ x) +
  facet_wrap(~Beach)


#' ## LME4 models 
#'
#' Random intercept model
lmer.ri <- lmer( ) #lme4 package
summary(lmer.ri)

#' Try to write down an equation describing the model. This will be easiest
#' to explore using an R Markdown (.Rmd) document open in a separate
#' window.


#' Random intercept and slope model
lmer.rc <- lmer( )
summary(lmer.rc)


#' Try to write down an equation describing the model.


#' Confidence intervals (see Section 18.12.4) for a comparison of methods
#' The default is method = "profile" (profile likelihood). Other options
#' are method = "boot" or method = "wald" (based on a Normal sampling
#' distribution).
confint(fit.mixed)


## -------------------------------------------------------------------
#'
#' ## Best-Linear Unbiased Predictions (BLUPs) for cluster-specific parameters 
#' 
#' Random effects = $b_i$ (deviations from mean parameter)
ranef(lmer.rc)
broom.mixed::tidy(lmer.rc, effects="ran_vals")


# fig.align='center', out.width="100%", fig.height=6, fig.width=12
plot_model(lmer.rc, type="re")


#' Cluster-specific parameters = fixed effects + random effects 
coef(lmer.rc)

#' For estimates of uncertainty, bootstrap is recommended (demonstrated
#' below for random slopes)
set.seed(101)
boot.pred <- bootMer(lmer.rc,
                     FUN=function(x){
                       rancoef<-coef(x)$Beach[,3]},
                     nsim=500)   


dim(boot.pred$t)
predboot.CI<-apply(boot.pred$t,2, quantile, c(0.025, 0.975))
predboot.CI


## -------------------------------------------------------------------------------------

#' 
#' ## Comparable Fixed effects models 
#' 
#' Beach specific intercepts
RIKZdat$Beach <- as.factor(RIKZdat$Beach)
lm.fc1 <- lm(Richness ~ Beach + NAPc - 1, data = RIKZdat)
summary(lm.fc1)

lm.fc1b <- lm(Richness ~ Beach + NAPc + exposure.c -1, data = RIKZdat)
summary(lm.fc1b)


##  Beach specific slopes
lm.fc2 <- lm( )
summary(lm.fc2)


## -------------------------------------------------------------------------------------

#' Compare fixed versus random effects 
#' First, let's refit the mixed effects model, but without exposure.c so that
#' the models are directly comparable
lmer.rc2 <- lmer(Richness ~ NAPc + (1 + NAPc | Beach), data = RIKZdat)
coef.fixed <- coef(lm.fc2)
coef.random <- coef(lmer.rc2)$Beach 
coefdat <- data.frame(intercepts = c(coef.fixed[1:9], coef.random[,1]), 
                      slopes = c(coef.fixed[10:18], coef.random[,2]), 
                      method = rep(c("Fixed", "Random"), each=9),
                      Beach = rep(names(coef.fixed)[1:9], 2))


## fig.align='center', fig.cap="Comparison of fixed versus random effects parameters demonstrating the shrinkage property of random-effects.", out.width = "100%", fig.height=6, fig.width=12----
p1 <- ggplot(coefdat, aes(Beach, intercepts, col=method)) + 
  geom_point(size=2)+ ggtitle("Intercepts")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_hline(yintercept=summary(lmer.rc2)$coef[1,1])
p2 <- ggplot(coefdat, aes(Beach, slopes, col=method)) + 
  geom_point(size=2)+ ggtitle("Slopes")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_hline(yintercept=summary(lmer.rc2)$coef[2,1])

ggpubr::ggarrange(p1, p2, ncol=2, common.legend=TRUE, legend="bottom")        


## fig.align='center', out.width="75%", fig.height = 4, fig.width=6----
ggplot(coefdat, aes(intercepts, slopes, group=Beach, col=method)) + 
  geom_point(size=2) + geom_path() +
  geom_hline(yintercept=summary(lmer.rc2)$coef[2,1], lty = 2) +
  geom_vline(xintercept=summary(lmer.rc2)$coef[1,1], lty = 2)

 
#' Repeat but for a model where the intercepts and slopes are not
#' assumed to covary
lmer.rc2b <- lmer(Richness ~ NAPc + (1  | Beach) + (0 + NAPc | Beach), 
                  data = RIKZdat)
coef.fixed <- coef(lm.fc2)
coef.random <- coef(lmer.rc2b)$Beach 
coefdat <- data.frame(intercepts = c(coef.fixed[1:9], coef.random[,1]), 
                      slopes = c(coef.fixed[10:18], coef.random[,2]), 
                      method = rep(c("Fixed", "Random"), each=9),
                      Beach = rep(names(coef.fixed)[1:9], 2))


## fig.align='center', fig.cap="Comparison of fixed versus random effects parameters demonstrating the shrinkage property of random-effects.", out.width = "100%", fig.height=6, fig.width=12----
p1 <- ggplot(coefdat, aes(Beach, intercepts, col=method)) + 
  geom_point(size=2)+ ggtitle("Intercepts")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_hline(yintercept=summary(lmer.rc2b)$coef[1,1])
p2 <- ggplot(coefdat, aes(Beach, slopes, col=method)) + 
  geom_point(size=2)+ ggtitle("Slopes")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_hline(yintercept=summary(lmer.rc2b)$coef[2,1])

ggpubr::ggarrange(p1, p2, ncol=2, common.legend=TRUE, legend="bottom")        


## fig.align='center', out.width="75%", fig.height = 4, fig.width=6----
ggplot(coefdat, aes(intercepts, slopes, group=Beach, col=method)) + 
  geom_point(size=2) + geom_path() +
  geom_hline(yintercept=summary(lmer.rc2b)$coef[2,1], lty = 2) +
  geom_vline(xintercept=summary(lmer.rc2b)$coef[1,1], lty = 2)


## ----------------------------------------------------------------

#' ## Compare Beach-specific regression lines
#'
#' By default, predict will give predictions that include the BLUPs, 
#' i.e., we get beach-level predictions
#+ fig.align='center', fig.height=8, fig.width=10, out.width="90%"
# Subject-specific regression lines from the mixed-effects model
# Could also use: predict(lmer.rc, re.form=~(1+NAPc|Beach))
RIKZdat$sspred.rc<-predict(lmer.rc)

# Subject-specific regression lines from the fixed-effects model
RIKZdat$fepred<-predict(lm.fc2)

# Plot
ggplot(RIKZdat, aes(NAPc, Richness)) +
  geom_point(size = 2) +
  geom_line(aes(NAPc, sspred.rc), lwd = 1.3, col = "red") +
  geom_line(aes(NAPc, fepred), lwd = 1.3, col = "blue") +
  facet_wrap( ~ Beach)


#' To get population average predictions, we need to use
#' re.form = ~ NA.
# fig.height=4, fig.width=6, out.width="90%" 
RIKZdat$papred <- predict(lmer.rc, re.form = ~0)

ggplot(RIKZdat, aes(NAPc, Richness, shape=exposure.c, col = Beach)) +
  geom_point(size=2)+
  geom_line(aes(NAPc, papred), lwd=1.3, col="red")


## ----------------------------------------------------------------

#' ## Diagnostics
#' 
#' Default plots are residual versus fitted values (incorporating 
#' the BLUPs)
#+ fig.align='center', out.width = "90%",  fig.height=4, fig.width=8
ri<-plot(lmer.ri, main = "Random intercept") # residual versus fitted
rc<-plot(lmer.rc, main = "Random intercept and slope")
grid.arrange(ri, rc, ncol=2)


#+ fig.align='center', out.width = "80%", fig.height=12, fig.width=8, warning=FALSE, message = FALSE----
check_model(lmer.ri)


#+fig.align='center', out.width = "80%", fig.height=12, fig.width=8, warning=FALSE, message = FALSE----
check_model(lmer.rc)

 