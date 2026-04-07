#' ---
#' title: "17-LinearMixedEffectsModels1.R"
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
knitr::opts_chunk$set(fig.align='center', out.width = "55%", fig.height = 5, fig.width =6)
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
data(Selake)
head(Selake)

#' Fit naive linear model assuming all observations are independent
#+ fig.align='center', out.width="100%", fig.height=4, fig.width=9 
fit.naive <- lm(Log_fish_Se ~ Log_Water_Se, data=Selake)
plot.naive <- ggplot(data=Selake, 
    aes(x=Log_Water_Se, y=Log_fish_Se, col=Lake))+
    ggtitle("Se in Fish vs Se Water with linear regression line")+
    xlab("log(Se) in the water\nPoints jittered")+ylab("log(Se) in fish")+
    geom_point(size=3, position=position_jitter(width=0.05))+ 
    geom_abline(intercept=coef(fit.naive)[1], slope=coef(fit.naive)[2])+
                    theme(legend.position="none") 
residplot <- ggplot(broom::augment_columns(fit.naive, Selake), 
                    aes(.fitted, .resid, col=Lake)) +
                    geom_point() +
                    geom_smooth() +
                    ggtitle("Residual vs. fitted values") +
                    xlab("Fitted values") +
                    ylab("Residuals") 
ggpubr::ggarrange(plot.naive, residplot, ncol=2, common.legend=TRUE, legend="bottom")


#' Fit a model relating mean log(fish_Se) to log(water_Se) 
fishse.avg <- Selake %>% group_by(Lake) %>%
   dplyr::summarize(Log_Water_Se=mean(Log_Water_Se), fish.avg.se=mean(Log_fish_Se))
fit.avg <- lm(fish.avg.se ~ Log_Water_Se, data=fishse.avg)

#' Fit a model where each lake has its own random intercept
fit.mixed <- lme4::lmer(Log_fish_Se ~ Log_Water_Se + (1|Lake), data=Selake)

#' Compare results
modelsummary(list("Naive linear model" = fit.naive, 
                  "Linear model fit to averages" = fit.avg, 
                  "Mixed-effects model" = fit.mixed),
             gof_omit = ".*", estimate  = "{estimate} ({std.error})", statistic=NULL,  
             coef_omit = "SD",
             title="Regression coefficients (SE) from a naive linear regression model 
             assuming independence, a regression model fit to lake averages, and a 
             random-intercept model.")


#' Summary of mixed effects models
summary(fit.mixed)
confint(fit.mixed)


#' Note, we cannot fit a model with both lake Se and Lake level intercepts
#' using fixed effects when we have only 1 observation per lake.  The Se effect
#' is confounded with the lake level effects.
lm.fixed <- lm(Log_fish_Se ~ Log_Water_Se + Lake -1, data=Selake)
summary(lm.fixed)

