#' ---
#' title: "Non-linear Regression"
#' author: ""
#' format: 
#'  html:
#'    toc: true
#'    embed-resources: true
#' ---
#' 
#' 
#' Load libraries
#' 
#+ warning=FALSE, message=FALSE
library(knitr) # for reproducible report
library(Data4Ecologists) # for data
library(ggplot2) # for plotting
library(splines)  # for splines
library(performance) # for residual plots
library(dplyr) # for data wrangling
library(car) # for added variable plots
library(ggeffects) # for effect plots
#' 
#' Settings for Knitr (optional)
opts_chunk$set(fig.width = 8, fig.height = 6, out.width = "60%", fig.align = 'center')

#' ## Exploratory data analysis
#' 
#' 1. Access the data and look at first 6 records.
data(clutch) 
head(clutch)

#' 
#' 2. Change names to lower case
#' 
names(clutch) <- tolower(names(clutch))

#' 
#' 3. Change year to a factor variable and then plot the data.
#' 
clutch$year <- as.factor(clutch$year)

#' Plot the data
#+ fig.alt = "Scatterplot of clutch size veruss nest initiation date for each year"
ggplot(clutch, aes(x = date, y = clutch, color = year)) + geom_point() 

#' 
#' Or, with facets:
#' 
#+  fig.height=4, fig.width = 12, out.width = "100%", fig.alt = "Scatterplot of clutch size veruss nest initiation date for each year using facets for year"
ggplot(clutch, aes(x = date, y = clutch)) + geom_point() + facet_wrap(~ year)

#'
#' 4. Drop outliers & plot again. Note: if we were going to use these data in other programs, we would want to save them in our data folder under a different name & document how the data differ from the raw data set.
#' 
clutch.r <- subset(clutch, clutch <= 15) #clutch.r for reduced

#' Plot again
#+ fig.height=4, fig.width = 12, out.width = "100%",  fig.alt = "Scatterplot of clutch size veruss nest initiation date for each year using facets for year"
ggplot(clutch.r, aes(x = date, y = clutch)) + geom_point() + facet_wrap(~ year)

#' 
#' 
#' 5. Ignore the fact that most of the nest structures were observed for
#' multiple years (i.e., the observations are not necessarily independent).  Fit a linear model relating `clutch` to `year` and `date`.  Create a residual versus fitted value plot for these data. Comment on whether you think the assumptions of linear regression are reasonably met.
#'  
lm.fit1 <- lm(clutch ~ year + date, data = clutch.r)

#' 
#' Let's use the `check_model` function in the performance package:
#' 
#+ fig.height=8, fig.width=8, out.width = "80%", warning=FALSE, message=FALSE,  fig.alt = "Residual plots"
check_model(lm.fit1, check = c("linearity", "homogeneity", "qq", "normality"))

#' Now, lets let at added variable and component + residual plots
#+ fig.height=8, fig.width=8, out.width = "80%", fig.alt = "Added variable and component + residual plots"
avPlots(lm.fit1, terms = "date")
crPlots(lm.fit1, terms = "date")

 
#' We see that there is a clear non-linear trend in the residuals versus predicted values plot (top left) suggesting that the linearity assumption is problematic. In addition, the residuals fall away from the qqline at the extremes, suggesting that the errors are not Normally distributed (the distribution has longer "tails" than the Normal distribution). You probably also notice the weird patterning in the residual plots.  This is due to clutch sizes taking on only integer values.  Because these are count data, we might want to consider a distribution other than a Normal distribution (e.g., using methods we will learn about later in the book).


#'
#' ## Polynomial model
#'  
#' 6. Fit a model using `lm`, allowing for a quadratic relationship between `date` and `clutch`. Again, create a residual versus fitted value plot.  Evaluate whether the assumptions of linear regression are reasonably met.
#'  
lm.fit2<-lm(clutch ~ year + poly(date, 2), data = clutch.r)

# Component + residual plot
#+  fig.alt = "Added variable and component + residual plots"
par(mfrow=c(1,2))
termplot(lm.fit2, partial.resid = TRUE, se = TRUE)
par(mfrow=c(1,1))

#' Let's look at residual plots again
#+ fig.height=8, fig.width=8, out.width = "80%", fig.alt = "Residual plots"
check_model(lm.fit2, check = c("linearity", "homogeneity", "qq", "normality"))

#' 
#' The residuals versus fitted value plot looks much better. The qqplot and histogram of residuals suggest the Normality assumption might be slightly suspect (as might be expected since we have count data, which are clearly not continuous).
#' 

#' ## Linear spline
#' 
#' 7. Now, fit a linear spline model, with a single knot at `date` = 150. will need to create a new predictor that is 0 if `date` $<=150$ and = `date`-150 otherwise.
#' 
clutch.r$date.sp <-0
clutch.r$date.sp[clutch.r$date > 150] <- clutch.r$date[clutch.r$date > 150] - 150

#' 
#' Some other options for accomplishing the same task include:
## ---------------------------------------------
clutch.r$date.sp2<-ifelse(clutch.r$date < 150, 0, clutch.r$date -150)
clutch.r<-clutch.r %>% mutate(date.sp3 = ifelse(date < 150, 0, date -150))
all.equal(clutch.r$date.sp, clutch.r$date.sp2)
all.equal(clutch.r$date.sp, clutch.r$date.sp3)

#' 
#' Fit the linear spline model and inspect residuals.  These also look good!
#' 
lm.fit3<-lm(clutch ~ year + date + date.sp, data = clutch.r)

#' 
#' 
#+ fig.height=8, fig.width=8, out.width = "100%", fig.alt = "Residual plots"
check_model(lm.fit3, check = c("linearity", "homogeneity", "qq", "normality"))

#' 
#' Everything looks pretty good here.  There is no significant trend in the residual versus fitted value as we move from left to right. There is one large outlier, but otherwise, the Normality assumption seems reasonable. 
#' 
#' 8. Create a plot showing the fit of the model, overlaid on the original data.  Lets get predictions for date = 90 to 180.
#' 
newdata <- data.frame(expand.grid(date = seq(90,180, by = 1), 
                                year = c("1997", "1998", "1999")))
newdata$date.sp <- ifelse(newdata$date <= 150, 0, newdata$date - 150)
newdata$fitted<-predict(lm.fit3, newdata = newdata)

#' 
#' Plot fit of model overlaid on the data using plots and lines:
#+ fig.alt = "plot of spline model overlayed on data"
ggplot(newdata, aes(date, fitted, color = year)) + geom_line() +
  geom_point(data = clutch.r, aes(date, clutch)) + xlab("Nest Initiation date") +
   ylab("Fitted values") + theme(axis.title = element_text(size = 16))


#' 
#' Or, with facets
#' 
#+ fig.height=4, fig.width = 12, out.width = "100%", fig.alt = "plot of spline model overlayed on data with facets for year"
ggplot(newdata, aes(date, fitted, color = year)) + geom_line()+ 
  geom_point(data = clutch.r,aes(date, clutch)) + 
  xlab("Nest Initiation date") + ylab("Fitted values") + facet_wrap(~ year) + 
  theme(axis.title = element_text(size=16))

#' 
#' ## Cubic spline model
#' 
#' 9. Fit a model with natural cubic regression splines and 3 degrees of freedom allocated to date.
#' 
#' Spline Fit and inspect residuals.  This also looks good.
#' 
lm.fit4<-lm(clutch~ year + ns(date, 3), data = clutch.r)

#+fig.height=8, fig.width=8, out.width = "80%", fig.alt = "Residual plots"
check_model(lm.fit4, check = c("linearity", "homogeneity", "qq", "normality"))

#' Predictions for each year, mean date using `ggpredict` package
mu.hat <- ggpredict(lm.fit4)
print(mu.hat)

#' Predictions weighted and averaged across years and for mean date
mu.hat.marg <- ggeffect(lm.fit4)
print(mu.hat.marg)

# Both have plotting methods that can be used to visualize the results
#+ fig.alt = "Effect plots for cubic spline model"
plot(mu.hat)
plot(mu.hat.marg)


#' Plot of fitted model.  lets get predictions for date = 90 to 180
newdat<-data.frame(expand.grid(date = seq(90,180, by=1), 
                               year = c("1997","1998","1999")))
newdat$fitted <- predict(lm.fit4, newdata = newdat)

#+ fig.height=4, fig.width=12, out.width="100%", fig.alt = "Plot of model predictions."
ggplot(newdat, aes(date, fitted)) + geom_line() + 
  geom_point(data = clutch.r,aes(date, clutch)) + 
  xlab("Nest Initiation date") + ylab("Clutch size") + facet_wrap(~ year) +
  theme(axis.title = element_text(size=16))

 
#' ## Footer
#' 
#' We can render the document using this code (you will need to modify it for your directory structure):
#' 
#' quarto::quarto_renderr(here("RCode", "04-1Non-LinearRegression.R"))
#' 
#' Session Information:
#' 
sessionInfo()
 
