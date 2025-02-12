#' # Multicollinearity 
#' 
#'  
#+ warning=FALSE, message=FALSE
library(dplyr) # for data wrangling

 
#' ## Motivating example:  what factors predict how long mammals sleep? 
#' 
#' Load sleep data, drop observations with missing variables and look at 
#' a pairwise scatterplot 
library(openintro)
data(mammals, package="openintro") 
mammals <- mammals %>% dplyr::select(total_sleep, life_span, gestation,
                          brain_wt, body_wt, predation, 
                          exposure, danger) %>%
  filter(complete.cases(.)) 
GGally::ggpairs(mammals)

#' Lets fit some different models and explore the estimates and SEs
library(modelsummary)
model1<-lm(total_sleep ~ life_span, data=mammals)
model2<-lm(total_sleep ~ life_span + danger + log(brain_wt), data=mammals)
model3 <- lm(total_sleep ~ life_span + gestation + log(brain_wt) + 
               log(body_wt) + predation + exposure + danger, data=mammals)
model4 <- lm(total_sleep ~ predation, data=mammals)
modelsummary(list("Model 1" = model1, "Model 2" = model2, "Model 3" = model3, "Model 4" = model4), gof_omit = ".*",
             estimate = "{estimate} [{std.error}] ({p.value})",
             statistic = NULL, 
             title="Estimates [SE] (p-values) for different models.")


mosaic::cor(life_span ~ log(brain_wt), data = mammals, use = "complete.obs")


#' ## Variance inflation factors (VIF)
car::vif(model3)

#'Visualization of VIFS
performance::check_model(model3, check = "vif")

