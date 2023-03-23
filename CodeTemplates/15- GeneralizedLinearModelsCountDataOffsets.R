#' ---
#' title: "Offsets"
#' author: "John Fieberg"
#' date: "2023-03-15"
#' output: 
#'    html_document:
#'      toc: true
#'      toc_depth: 2
#'      toc_float: true
#' ---

#' Load data
#+ warning=FALSE, message=FALSE
library(Data4Ecologists)
data(reeffish)
head(reeffish)


#' Create log(Duration) variable to measure survey effort 
reeffish <- reeffish %>% mutate(logduration = log(Duration))
fitoffset <- glm.nb(Seabass ~ as.factor(Year) + poly(Temperature,2) + 
                      Depth + Doy + Hour + offset(logduration) -1, data = reeffish)
summary(fitoffset)


#' Fit the same model but estimate the coefficient for logduration
fitpredictor <- glm.nb(Seabass ~ as.factor(Year) + poly(Temperature,2) + 
                         Depth + Doy + Hour + logduration -1, data = reeffish)
summary(fitpredictor)

#' Profile-likelihood confidence interval
round(confint(fitpredictor, parm = "logduration"), 2)
AIC(fitoffset, fitpredictor)
  