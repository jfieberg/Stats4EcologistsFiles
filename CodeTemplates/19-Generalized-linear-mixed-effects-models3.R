#' ---
#' title: "19-Generalized-linear-mixed-effects3.R"
#' author: "John Fieberg"
#' output: 
#'    html_document:
#'      toc: true
#'      toc_depth: 2
#'      toc_float: true
#' ---
#' 

#' Load lme4 package
#+ warning=FALSE, message = FALSE
library(lme4)
library(dplyr)
library(ggplot2)
library(gridExtra) # for multi-panel plots

#' This program considers the simulation example in the book (CH 19) 
#' The random intercept, $b_{0i}$, accounts for among-individual differences
#' in their initial abilities and experiences with skating. 
#' Let's simulate 10 observations for each of 100 individuals, setting:
#' 
#' - $\beta_0 = -6$
#' - $\beta_1 = 1$
#' - $\sigma_{b_0} = 2$
#' 
## --------------------------------------------------------
# ## Simulation set up
set.seed(1200)

# Sample size
nindividuals <- 100 # number of individuals
nperindividual <- 10 # number of observations per individual
ntotal <- nindividuals * nperindividual # total number of observations

# Data generating parameters
b0_i <- rnorm(nindividuals, 0, 2) # deviations from mean intercept
beta_0 <- -6 # mean intercept
beta_1 <- 1 # slope for number of practices

# Generate 10 different observations for each individual
# for practices 1 through 10
practice.num <- seq(1, 10, length = 10)
individual <- rep(1:nindividuals, each = nperindividual) # individual id

# Y|boi is Bernouli p[ij]
logitp <- beta_0 + rep(b0_i, each = nperindividual) + practice.num * beta_1
p <- plogis(logitp)
Y <- rbinom(ntotal, size = 1, prob = p)

# Create data set for model fitting
skatedat <- data.frame(Y = Y, practice.num = practice.num, individual = individual)

#' 
#' We can then fit the GLMM used to generate the data and show that we can recover the simulation parameters:

glmermod <- glmer(Y ~ practice.num + (1 | individual), family = binomial(),
                  data = skatedat)
summary(glmermod)

 
#' For example, the true parameters are contained within Wald-based confidence
#' intervals formed using a Normality assumption.
confint(glmermod, method="Wald")

#' 
#' Now, let's plot the average value of $Y$ at each time point - i.e., the 
#' proportion of the 100 skaters that can successfully skate backward at each 
#' practice.  We will also plot the predicted mean response curve formed by 
#' setting the random effect, $b_{0i} = 0$. The latter response curve 
#' corresponds to the predicted probability that a "typical skater" (i.e.,
#'  the median individual) will be able to skate backwards at each practice:
#' 
#+ fig.align='center', fig.align='center', fig.width=6, fig.height=4, out.width="75%"
# Predictions for a "typical subject" with b0i=0
practice.num <- seq(1, 10, length = 10)
mod <- data.frame(practice.num = practice.num)
mod$glmer.logit <- fixef(glmermod)[1] + practice.num * fixef(glmermod)[2]
mod$glmer.p <- plogis(mod$glmer.logit)

# Summarize average Y (i.e., proportion of skaters that successfully
# skate backwards) at each practice 
mdat <- skatedat %>%
  group_by(practice.num) %>%
  dplyr::summarize(PopProportion = mean(Y))

# Combine two curves and plot 
mdat2 <- data.frame(
  practice.num = rep(practice.num, 2),
  p = c(mod$glmer.p, mdat$PopProportion),
  type = rep(c("Probability, typical individual", "Sample proportion"), each = 10)
)

# Plot
ggplot(mdat2, aes(practice.num, p, col=type)) +
  geom_line() + geom_point() +
  xlab("Number of practices") +
  ylab("Proportion of skaters that can skate backwards")


#' We see that the two curves do not line up - the predicted probability of skating backwards for a "typical individual" (with $b_{0i} = 0$) differs from the proportion of the sample that can skate backwards. This occurs because of the non-linear logit transformation, $E[Y | X] \ne  E[Y | X, b_{0i} = 0]$.  To understand this mismatch, we next compare:
#' 
#' - the mean response curve for each individual on both the logit and probability scales 
#' - the mean of these individual-specific curves on both the logit and probability scales
#' - the overall mean population response curve on both the logit and probability scales
#' 
#' We begin by computing the individual-specific curves, below:
#' 
## --------------------------------------------------------

pdat <- NULL
for (i in 1:nindividuals) {
  logitp.indiv <- beta_0 + b0_i[i] + practice.num * beta_1 # individual i's curve on logit scale
  p.indiv <- plogis(logitp.indiv) # individual i's curve on the probability scale
  tempdata <- data.frame(p = p.indiv,
                         logitp = logitp.indiv,
                         individual = i)
  pdat <- rbind(pdat, tempdata)
}   
pdat$practice.num<-practice.num

#' 
#' We then calculate the average of the individual-specific curves on both the logit and probability scales:
#' 
#' 
## --------------------------------------------------------
pop.patterns<-pdat %>% group_by(practice.num) %>% 
  dplyr::summarize(meanlogitp=mean(logitp), meanp=mean(p))

#' 
#' 
#' Lastly, we compare the individual-specific (i.e., subject-specific or 
#' conditional means; in black) to their average (in red) and to the 
#' population-average response curves (in blue) on both logit and probability 
#' response scales:
#' 
#+ fig.align='center', fig.width=8, fig.height=4, out.width="95%"
m1 <- ggplot(pdat, aes(practice.num, logitp)) +
  geom_line(aes(group = individual)) + ylab("logit p") +
  geom_line(data = pop.patterns, aes(practice.num, meanlogitp), col = "blue", size = 1.2) +
  geom_line(data = mod, aes(practice.num, glmer.logit), col = "red", size = 1.2, linetype = "dashed") +
  theme(legend.position = 'none') + xlab("Practice Number") +
  ylab(expression(Logit(p[i]))) + ggtitle("A)")

m2 <- ggplot(pdat, aes(practice.num, p)) +
  geom_line(aes(group = individual)) + ylab("p") +
  geom_line(data = pop.patterns, aes(practice.num, meanp), col = "blue", size = 1.2) +
  geom_line(data = mod, aes(practice.num, glmer.p), col = "red", size = 1.2, linetype = "dashed") +
  theme(legend.position = 'none') + xlab("Practice Number") +
  ylab(expression(p[i])) + ggtitle("B)")
grid.arrange(m1, m2, ncol = 2)


#' To calculate the population mean, we need to integrate over the random effects.
#+fig.align='center', fig.width=6, fig.height=4, out.width="75%"
# Matrix to hold marginal mean response curve, E[Y|X] 
pa.rate<-matrix(NA,10,1)  
sigma2b0 <- as.data.frame(VarCorr(glmermod))[1, 4]
for(i in 1:10){
  intfun<-function(b0){
    plogis(fixef(glmermod)[1] + b0 + practice.num[i]*fixef(glmermod)[2])*dnorm(b0,0,sqrt(sigma2b0))}
  pa.rate[i]<-integrate(intfun,-Inf, Inf)[1]
}              
mdat4<-rbind(mdat2, 
             data.frame(practice.num = practice.num,
                        p = unlist(pa.rate),
                        type = rep("Marginal mean by Integration", 10)))

# Plot
ggplot(mdat4, aes(practice.num, p, col=type)) +
  geom_line() + geom_point() +
  xlab("Number of practices") +
  ylab("Proportion of skaters that can skate backwards")

#' 
#' ### Approximating marginal means using GLMMadpative
#'
#+ warning = FALSE, message=FALSE
library(GLMMadaptive)
fit.glmm <- mixed_model(fixed = Y ~ practice.num,
                        random = ~ 1 | individual,
                        family = binomial(link = "logit"),
                        data = skatedat
)
summary(fit.glmm)
summary(glmermod)

marginal_beta <- marginal_coefs(fit.glmm, std_errors = FALSE)
mdat5<-rbind(mdat4, 
             data.frame(practice.num = practice.num,
                        p = plogis(marginal_beta$betas[1]+practice.num * marginal_beta$betas[2]),
                        type = rep("GLMMadaptive", 10)))

# Plot
ggplot(mdat5, aes(practice.num, p, col=type)) +
  geom_line() + geom_point() +
  xlab("Number of practices") +
  ylab("Proportion of skaters that can skate backwards")


#' Marginal odds ratios: the odds of skating backwards increases by a factor of 
#' 1.8 at the population level for each additional practice. 
exp(marginal_beta[1]$betas[2])

#' Compare to subject specific odds ratio (which will always be bigger) 
exp(fixef(glmermod)[2])

 