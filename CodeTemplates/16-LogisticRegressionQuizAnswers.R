#' Answers to google quiz

#' Read in a moose data set  
library(SightabilityModel)
data(exp.m)
str(exp.m)

#' Fit mod2
exp.m$year <- as.factor(exp.m$year)
mod2 <- glm(observed ~ voc + year, data = exp.m, family = binomial())
summary(mod2)

#' We could calculate things "by hand"
#' 2007: log odds
lp1<-coef(mod2)[1]+coef(mod2)[2]*0+coef(mod2)[4]

# 2005:  log odds
lp2<-coef(mod2)[1]+coef(mod2)[2]*0

# Then, calculate p using either
plogis(lp1); # or
exp(lp1)/(1+exp(lp1))


#' Or, we could use the predict function
newdat<-data.frame(voc=c(10, 0, 0, 50, 50), 
                   year = c("2006", "2005", "2007", "2005", "2007"))

#' Predictions for p
phats<-predict(mod2, newdata = newdat, type="resp")                   

#' Predictions for logit p
logitphats<-predict(mod2, newdata = newdat, type="link")

#Answers
phats[1]
logitphats[3]-logitphats[2]
logitphats[5]-logitphats[4]

phats[3]-phats[2]
phats[5]-phats[4]
