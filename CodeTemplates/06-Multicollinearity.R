#' # Multicollinearity 
#' 
#' 
#' **Learning objectives**
#' 
#' 1. Be able to describe and identify different forms of multicollinearity.
#' 2. Understand the effects of multicollinearity on parameter estimates and their standard errors. 
#' 3. Be able to evaluate predictors for multicollinearity using variance inflation factors.
#' 4. Understand some strategies for dealing with multicollinearity.
#' 
#' ## R Packages
#' 
#' We begin by loading the `dplyr` package:
#' 
## ----warning=FALSE, message=FALSE---------------------
library(dplyr) # for data wrangling

 
#' ## Motivating example:  what factors predict how long mammals sleep? 
#' 
#' Load sleep data, drop observations with missing variables and look at 
#' a pairwise scatterplot 
## ----sleeppairs, comment=NA, warning=FALSE, message=FALSE, fig.height=8, fig.width = 8, fig.align='center', fig.cap = "(ref:gallysleep)", out.width = "95%"----
library(openintro)
data(mammals, package="openintro") 
mammals <- mammals %>% dplyr::select(total_sleep, life_span, gestation,
                          brain_wt, body_wt, predation, 
                          exposure, danger) %>%
  filter(complete.cases(.)) 
GGally::ggpairs(mammals)

## ----sleep1, comment=NA, echo=c(-1), warning=FALSE, message=FALSE----
library(modelsummary)
model1<-lm(total_sleep ~ life_span, data=mammals)
model2<-lm(total_sleep ~ life_span + danger + log(brain_wt), data=mammals)


## ----sleep2-------------------------------------------
model3 <- lm(total_sleep ~ life_span + gestation + log(brain_wt) + 
               log(body_wt) + predation + exposure + danger, data=mammals)
model4 <- lm(total_sleep ~ predation, data=mammals)
modelsummary(list("Model 1" = model1, "Model 2" = model2, "Model 3" = model3, "Model 4" = model4), gof_omit = ".*",
             estimate = "{estimate} [{std.error}] ({p.value})",
             statistic = NULL, 
             title="Estimates [SE] (p-values) for different models.")


## ----comment=NA, cache.lazy=FALSE, warning=FALSE, message=FALSE----
mosaic::cor(life_span ~ log(brain_wt), data = mammals, use = "complete.obs")


#' ## Variance inflation factors (VIF)
car::vif(model3)

#' 
#' We can also use the `check_model` function in the performance package to create a nice visualization of the VIFS (Figure \@ref(fig:vifcheck)). We see that several of the VIFs are large and greater than 10, suggesting that several of our predictor variables are collinear.  You will have a chance to further explore this data set as part of an exercise associated with this Section.
#' 
#' (ref:checkvif) Variance inflation factors visualized using the `check_model` function in the `performance` package [@performance].
## ----vifcheck, fig.cap = "(ref:checkvif)", fig.align='center', out.width="100%", fig.height=6, fig.width=12----
performance::check_model(model3, check = "vif")

