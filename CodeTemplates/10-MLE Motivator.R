#' Motivation for maximum likelihood
#' 
#' #+warning=FALSE, message=FALSE
library(Data4Ecologists)
library(ggplot2)
library(performance)

#' Read in data  
data(slugs)

#' Test for differences in density among the two sites using a t-test
t.test(slugs~field, data=slugs)


#' Assumes data are normally distributed within each group
boxplot(slugs~field, data = slugs)

#' Also, what about data sets where we have:
#' 
#' - multiple predictors
#' - counts (is a Normall distribution appropriate)
#' - variance will often increase with the mean


#' Longnose dace data
data("longnosedace")

#' Linear regression assumptions are not met!
#+ fig.width= 8, fig.height=8, out.width=80%, fig.alt="Residual plots"
lmdace <- lm(longnosedace ~ acreage + do2 + maxdepth + no3 + so4 + temp,
             data = longnosedace)
check_model(lmdace, check = c("linearity", "homogeneity", "qq", "normality"))

#' Given that we have count data, one could try a glm 
#' 
slug.glm <- glm(slugs~field, data= slugs, family = poisson())
summary(slug.glm)

