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
library(lme4)  # for fitting random-effects models
library(glmmTMB) # for fitting random-effects models
library(ggplot2)# for plotting
library(performance) # for model diagnostics


#' Load data sets
#+warning=FALSE, message=FALSE 
library(Data4Ecologists) # for data
data("RIKZdatdat")

#' Only 1 beach has the lowest exposure level:  Recode exposure to have only 2 levels.	
RIKZdat$exposure.c<-"High"	  
RIKZdat$exposure.c[RIKZdat$exposure%in%c(8,10)]<-"Low"	
RIKZdat$NAPc = RIKZdat$NAP-mean(RIKZdat$NAP) #center NAP variable	


#' Fit Poisson-Normal models using glmer to RIKZdat.
#' First, fit the random intercept model
glmer.ri <- glmer( ) #lme4 package
summary(lmer.ri)
 
#' Now, fit the random intercept and slope model
glmer.rc <- glmer( ) #lme4 package
summary(glmer.rc)

#' Check assumptions using the performance package
check_model(glmer.ri)
check_model(glmer.rc)

#' Fit negative binomial models using glmmTMB
glmmtmb.ri <- glmmTMB(  , family = "nbinom2")
summary(glmmtmb.ri)

glmmtmb.rc <- glmmTMB(  , family = "nbinom2") 
summary(glmmtmb.ri)

#' Check assumptions
check_model(glmer.ri)
check_model(glmer.rc)

#' Compare using AIC
AIC(glmer.ri, glmer.rc, glmmtmb.ri, glmmtmb.rc)
