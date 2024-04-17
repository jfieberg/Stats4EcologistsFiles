#' ---
#' title: "17-LinearMixedEffectsModels3.R"
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

#' Only 1 beach has the lowest exposure level:  Recode exposure to have only 2 levels.	
RIKZdat$exposure.c<-"High"	  
RIKZdat$exposure.c[RIKZdat$exposure%in%c(8,10)]<-"Low"	
RIKZdat$NAPc = RIKZdat$NAP-mean(RIKZdat$NAP) #center NAP variable	

#' Models fit using lme 
lme.fit<-lme(Richness~NAPc+exposure.c, random=~1|Beach, data=RIKZdat)
summary(lme.fit)

 
#' Degrees of freedom approximation using lme
nrow(RIKZdat) - length(unique(RIKZdat$Beach))- 1
length(unique(RIKZdat$Beach))- 1 -1


#' Better options are available (and implemented) when lmerTest has been
#' lodaed.
#+ warning=FALSE, message=FALSE
library(lmerTest)
lmer.ri <- lmer(Richness ~ exposure.c + NAPc + (1|Beach), data=RIKZdat) #lme4 package
summary(lmer.ri)
summary(lmer.ri,  ddf = "Kenward-Roger")


#' Loading lmerTest also changes the way `anova` works (it uses different
#' degrees of freedom calculations and no longer implements sequential tests)
anova(lmer.ri)

#' Compare results for random intercept versus random slope model
#' for predictors that do and do not vary withing a beach (exposure.c versus 
#' NAPc) 
lmer.rc <- lmer(Richness ~ exposure.c + NAPc + (1+NAPc | Beach), data=RIKZdat)
summary(lmer.rc)

#' Simulation-based likelihood ratio test
lrsimtest <- pbkrtest::PBmodcomp(lmer.rc, lmer.ri, nsim = 500)
summary(lrsimtest)


#' AIC is possible, but several caveats apply - see Bolker's GLMM FAQ
AIC(lmer.ri, lmer.rc)


#' Jack Weiss's strategy 
## -------------------------------------------------------------------------------------
pooled.model <- lm(Richness ~ NAPc, data = RIKZdat)
uncond.means.model <- lmer(Richness ~ 1 + (1 | Beach), data = RIKZdat)
randint.model <- lmer(Richness ~ 1 + NAPc + (1 | Beach), data = RIKZdat)
randcoef.model <- lmer(Richness ~ 1 + NAPc + (NAPc | Beach), data = RIKZdat)
AIC(pooled.model, uncond.means.model, randint.model, randcoef.model)

#' So, random coefficient model, then consider exposure.c, but fit 
#' both models using ML rather than REML
lmer.ri.MLa <- lmer(Richness ~ 1 + NAPc + (NAPc | Beach), 
                 REML = FALSE, data = RIKZdat)
lmer.ri.MLb <- lmer(Richness ~ 1 + NAPc + exposure.c + (NAPc | Beach),
                 REML = FALSE, data = RIKZdat)
AIC(lmer.ri.MLa, lmer.ri.MLb)

#' What if we had used REML?
lmer.ri.REMLa <- lmer(Richness ~ 1 + NAPc + (NAPc | Beach), 
                 data = RIKZdat)
lmer.ri.REMLb <- lmer(Richness ~ 1 + NAPc + exposure.c + (NAPc | Beach),
                 data = RIKZdat)
AIC(lmer.ri.REMLa, lmer.ri.REMLb)

#+ warning=FALSE, message = FALSE
library(lmerTest)
anova(lmer.ri.REMLb) 

## ------------------------------------------------------------
#'
#' ## Marginal model
#' 
#' Marginal model fit using GLS 
gls.fit<-gls(Richness ~ NAPc + exposure.c, method="REML",
             correlation=corCompSymm(form= ~ 1 | Beach),
             data=RIKZdat)
summary(gls.fit)
fixef(lmer.ri)


## -------------------------------------------------------------------------------------
# Extract variance parameters from the random intercept model
(variancepars<-as.data.frame(VarCorr(lmer.ri)))

# rho calculated from random intercept model
variancepars[1,4]/(variancepars[1,4]+variancepars[2,4])

# residual standard error in marginal model
sqrt(variancepars[1,4]+variancepars[2,4])


## -------------------------------------------------------------------------------------
#'
#' ## Multiple crossed random effects
library(Data4Ecologists)
data(HRData)


# fig.align='center', out.width = "70%", fig.height = 8, fig.width =8
HRsummary <- HRData %>% group_by(PackID, Year) %>%
  count()
ggplot(HRsummary, aes(Year, PackID, size = n)) + geom_point() + 
  ylab("Pack ID")


#' Models with multiple random effects
model1 <- lmer(log(HRsize) ~ Season + StudyArea + DiffDTScaled + LFD*EVIScaled + 
                 (1 | AnimalId) + (1 | Year),  REML=TRUE, data = HRData)
summary(model1)

model2 <- lmer(log(HRsize) ~ Season + StudyArea + DiffDTScaled + LFD*EVIScaled + 
                 (1 | PackID) + (1 | Year),  REML=TRUE, data = HRData)
summary(model2)

model3 <- lmer(log(HRsize) ~ Season + StudyArea + DiffDTScaled + LFD*EVIScaled + 
                 (1 | AnimalId) + (1 | PackID ) + (1 | Year), REML=TRUE, data = HRData)

summary(model3)
AIC(model1, model2, model3)

# Model with "PackYear" 
HRData <- HRData %>% mutate(PackYear = paste(HRData$PackID, HRData$Year, sep=""))

model5 <- lmer(log(HRsize) ~ Season + StudyArea + DiffDTScaled + LFD*EVIScaled + 
                 (1 | PackYear), REML=TRUE, data = HRData)
summary(model5)
AIC(model2, model5)
 

 