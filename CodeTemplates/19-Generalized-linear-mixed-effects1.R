#' ---
#' title: "19-Generalized-linear-mixed-effects1.R"
#' author: "John Fieberg"
#' output: 
#'    html_document:
#'      toc: true
#'      toc_depth: 2
#'      toc_float: true
#' ---
#' 

#' Load libraries
#+warning=FALSE, message=FALSE
library(ggplot2);theme_set(theme_bw())
knitr::opts_chunk$set(fig.align='center', 
                      out.width = "55%", 
                      fig.height = 5, 
                      fig.width =6)
options(width=160)

library(tidyverse) # for data wrangling
library(gridExtra) # for multi-panel plots
library(lme4)  # for fitting random-effects models
library(nlme) # for fitting random-effect and gls models
library(glmmTMB) # for fitting random-effects models
library(modelsummary) # for creating output tables
library(kableExtra) # for tables 
options(kableExtra.html.bsTable = T)
library(ggplot2)# for plotting
library(performance) # for model diagnostics



#' Load clutch data
library(Data4Ecologists) # data package
data(clutch)

#' Get rid of observations with nest initiation dates > 30 May
#' and a few outliers representing nests that were likely parasitized
#+ fig.height=4, fig.width=12, out.width = "90%"
clutch <- clutch %>% filter(CLUTCH < 13 & date < 149) %>%
  mutate(year = as.factor(year)) %>% 
  mutate(Ideply = ifelse(DEPLY==2, TRUE, FALSE))
ggplot(clutch, aes(date, CLUTCH, col=Ideply))+
  geom_point() + 
  geom_smooth(method= lm, formula = y ~ x) +
  facet_wrap(~year)

#' In this data set, we have 2 potential grouping variables, `year` 
#' and `strtno`, the latter is a unique identifier for each structure.  
#' 
#' **Think-pair-share**: is it better to model the effect of these 
#' grouping variables (and any interactions with them) as fixed or 
#' random effects?  What is the difference?

#' Fit random intercept model
clutch.ri <- lmer(CLUTCH ~ date + Ideply + (1 | strtno), data=clutch)
summary(clutch.ri)

#' Fit marginal model using gls
clutch.gls <- gls(CLUTCH ~ date + Ideply,  
                        corCompSymm(value = 0.3, form = ~ 1 | strtno), data=clutch)

#' Compare estimates, standard errors, etc
modelsummary(list("Random intercept" = clutch.ri, "GLS" = clutch.gls),
             gof_omit = 'BIC|REMLcrit')  

#' Not sure what is going on with the AIC for the GLS model, but it is likely that 
#' a constant in the log-likelihood function is dropped in either gls or lmer. 
#' We get identical AIC's if we use the AIC function:
AIC(clutch.ri, clutch.gls)


#' Extract variance parameters from the random intercept 
VarCorr(clutch.ri)
variancepars<-as.data.frame(VarCorr(clutch.ri))
variancepars

#' Compare the intraclass correlation coefficient, $\rho$ to variance parameters
#' in the mixed model.
#' 
#' rho = $\sigma^2_b/(\sigma^2_b + \sigma^2_{\epsilon})$ 
variancepars[1,4]/(variancepars[1,4]+variancepars[2,4])
coef(clutch.gls$modelStruct$corStruct, unconstrained = FALSE)

#' Residual standard error in marginal model:
#' 
#' $\sigma^2 = \sqrt{\sigma^2_b + \sigma^2_{\epsilon}}$
sqrt(variancepars[1,4]+variancepars[2,4])
clutch.gls$sigma

 

#' ## Extensions to Count and Binary Data {#nomarg}
#' 
# load data from Data4Ecologists package
data(nestocc)

# turn deply, period, year, strtno into factor variables
nestocc$deply <- as.factor(nestocc$deply)
nestocc$period <- as.factor(nestocc$period)
nestocc$year <- as.factor(nestocc$year)
nestocc$strtno <- as.factor(nestocc$strtno)

#' *Fixed or random effects*: note we have 3 categorical predictors here:
#'   `period`, `year` and `strno`. Which should we include as random
#'   effects? 

# remove observations with missing records and fit model
nestocc <- nestocc %>% filter(complete.cases(.)) 
nest.ri <- glmer(nest ~ year + period + deply + vom +
      year * period + vom * period + poly(wetsize, 2) + (1 | strtno),
      family = binomial(), data = nestocc, control = glmerControl(optimizer = "bobyqa"))

 
#' Design matrix
nestocc[1:3,]
model.matrix(nest.ri)[1:3,]
summary(nest.ri)
car::Anova(nest.ri)

#' ## Marginal response patterns
#' 

#' #### Subject-specific (SS) and population-averaged (PA) relationships:
#' 
#'  Obtain model-based estimates of SS and PA occupancy rates at
#'  different levels of VOM for:
#'  
#'  - Year = 1997 (year1998=year1999=0)
#'  - period = 1 (period2-period4=0)
#'  - cylinder type = single (deply2=0)  
#'  - wetland size = 0 (and size^{2} = 0)
#'  - all interactions =0.
#'  
#'  Do this for the typical individual (random strtno = 0) and for individuals 1 and 2 standard deviations from the "typical" individual.
#'
#' Store fixed effect coefficients, sigma^, and set up matrix of predictors
betas <- fixef(nest.ri)
sigmab<- as.data.frame(VarCorr(nest.ri))[1,5]
vom<-seq(0.05, 2.28, length=100)  #VOM values 
Xcovs<-cbind(rep(1,100),rep(0,100),rep(0,100),rep(0,100),rep(0,100),rep(0,100),rep(0,100),vom, rep(0,100), 
             rep(0,100),rep(0,100),rep(0,100),rep(0,100),rep(0,100),rep(0,100),rep(0,100),rep(0,100),rep(0,100),
             rep(0,100))
plogitmean<-Xcovs%*%betas  # logit(phat)
taus<-c(-2,-1,0,1,2)*sigmab # random effects for SS means in Fig. 1

#' Determine SS relationships on logit scale
pinds<-matrix(rep(plogitmean,5),5,100, byrow=T)+taus 

#' Determine PA relationships using numerical integration.  Need to perform a separate integration for each unique value of VOM.
#' I.E.: for each value of VOM, to get marginal prediction, we have to integrate over the estimated
#'      distribution of the random effects (N[0,sigmab]).  
pa.rate<-matrix(NA,100,1)  # Matrix to hold PA relationship (expected occupancy rate as a function of VOM)
for(i in 1:100){
  intfun<-function(x){plogis(plogitmean[i]+x)*dnorm(x,0,sigmab)}
  pa.rate[i]<-integrate(intfun,-Inf, Inf)[1]
}

#' Plot marginal and conditional means as a function of VOM.
#+fig.width=8, fig.height=4   
#   windows(width=24, height=8) # open plotting window of reasonable size
par(mfrow=c(1,1),oma=c(1,1.5,1.5,1.5), mar=c(4.9, 4.1, 2.1, 3.1), oma=c(0,1,0,0),ask=F, cex.lab=1.4, cex.axis=1.2, 
    cex=1.2,cex.main=1.4, tcl=0.5, bty="l")  

# Plot on original scale, marginal relationship along with SS relationship
plot(vom,unlist(pa.rate), ylim=c(0,1), lwd=3, type="l",xlab="VOM", ylab="Probability of use", xlim=c(0,2), col="red")

# Now, SS relationships for subjects with random effects = 0, +/-1, +/- 2 standard deviations 
for(i in 1:5){
  if(i==3){lines(vom, plogis(pinds[i,]),lwd=2,col="blue", lty=2)} # Average or typical subject with re=0
  else{lines(vom, plogis(pinds[i,]),lwd=1, lty=2)}
}


# Illustrate Zeger et al.'s (1988) approximation to the marginal relationship 
betam<-(sigmab^2*0.346 +1)^(-1/2)*betas
plogitMean<-Xcovs%*%betam
lines(vom,plogis(plogitMean), col=alpha(gray(0), 0.9), lwd=3, lty=2)
 



 