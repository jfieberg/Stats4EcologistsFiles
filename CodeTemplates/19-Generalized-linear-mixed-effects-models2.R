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
data("RIKZdat")

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



#' ## Generalized Estimating Equations

#'
#'#' Load libraries (geepack) for model fitting, MuMIn for QIC function
#+warning=FALSE, message=FALSE  
library(geepack)
#library(MuMIn)

#' Full independence
glmind <- glm(Richness~ NAPc + exposure.c, family=poisson(), data=RIKZdat)
  
#' Independent working correlation 
gee.ind<-geeglm(Richness ~ NAPc + exposure.c, id=Beach, corstr="independence",
                family=poisson(), data=RIKZdat)
summary(gee.ind)

#' Exchangeable working correlation structure  
gee.exc<-geeglm(Richness ~ NAPc + exposure.c,id=Beach, corstr="exchangeable",
                family=poisson(), data=RIKZdat)
summary(gee.exc)

#' Compare fit  
QIC(gee.ind, gee.exc)

# Compare to glmer.rc  (coefficients )
summary(glmer.rc) 

#' Compare using ggplot (after constructing intervals for each approach)
ests<-c(coef(glmind), fixef(glmer.ri), fixef(glmer.rc), coef(gee.ind), coef(gee.exc))
se.gee.ind<-summary(gee.ind)$coef[,2]
se.gee.exc<-summary(gee.exc)$coef[,2]
cis.glm<-confint(glmind)
cis.glmer.rc<-confint(glmer.rc)
cis.glmer.ri<-confint(glmer.ri)
cis.glmer<-cis.glmer.rc[4:6,]
cis.glmer2<-cis.glmer.ri[2:4,]
lcl<-c(cis.glm[,1], cis.glmer2[,1], cis.glmer[,1], coef(gee.ind)-1.96*se.gee.ind, coef(gee.exc)-1.96*se.gee.exc)
ucl<-c(cis.glm[,2], cis.glmer2[,2], cis.glmer[,2], coef(gee.ind)+1.96*se.gee.ind, coef(gee.exc)+1.96*se.gee.exc)
parmest<-data.frame(ests=ests, ucl=ucl,lcl=lcl, parm=rep(names(gee.ind$coef),5), 
                    type=rep(c("glm","glmer.ri", "glmer.rc", "GEE.ind", "GEE.exc"), each=3))

ggplot(parmest, aes(x=type, y=ests, group=type)) +
  geom_errorbar(width=.1, aes(ymin=lcl, ymax=ucl), colour="blue") + 
  geom_point(shape=21, size=3, fill="white")+
  facet_wrap(~parm, ncol=2, scales="free")+ylab("Estimate")
