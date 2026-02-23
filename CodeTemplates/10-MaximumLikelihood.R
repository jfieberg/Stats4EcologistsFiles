#' # Maximum likelihood {#mle}
#' 
## ----echo=FALSE-------------------------------------------------
knitr::opts_chunk$set(cache=1, 
                      cache.extra = knitr::rand_seed, 
                      autodep = TRUE)
knitr::opts_chunk$set(out.width = "55%", 
                      fig.height = 5, 
                      fig.width = 6, 
                      fig.align = "center") 


#' ## R packages
#' 
#' We begin by loading a few packages upfront:
#' 
## ----warning=FALSE, message=FALSE-------------------------------
library(kableExtra) # for creating tables 
options(kableExtra.html.bsTable = T)
library(dplyr) # for data wrangling
library(ggplot2) # for plotting
library(patchwork) # for creating multi-panel plots

 
#' In addition, we will use data and functions from
#' the following packages:
#' 
#' - `Data4Ecologists` for `beargrowth`  
#' - `emdbook` [@emdbook] for the `deltavar` function used to 
#' calculate standard errors via the delta method.
## ----bearplot, warning=FALSE, message=FALSE, fig.cap="Weight versus age depicted for black bears monitored in Minnesota."----
library(Data4Ecologists)
theme_set(theme_bw())
data(beargrowth)
beargrowth <- beargrowth %>% filter(sex=="Female")
ggplot(beargrowth, aes(age, wtkg)) + 
  geom_point() + geom_smooth() + ylab("Weight (kg)") + 
  xlab("Age (years)") 

#' We can describe the model we want to fit using the following set of
#' equations:  
#' \begin{equation}
#' \begin{split}
#' Weight_i \sim N(\mu_i, \sigma^2_i)\\
#' \mu_i = L_{\infty}(1-e^{-K(Age_i-t_0)})\\
#' \sigma^2_i = \sigma^2 \mu_i^{2\theta}
#' \end{split}
#' \end{equation}
#' 

#' Likelihood function
## ---------------------------------------------------------------
logL<-function(pars, dat){

    mu<-linf*(1 - exp(-k*(dat$age - t0))) 
    vari<-sig^2*abs(mu)^(2*theta)
    sigi<-sqrt(vari)
    ll<- # Log likelihood
    return(ll) 
  }   

#' 
#' To get starting values, note:
#' 
#' - The smooth curve asymptotes at around 100, so we use
#'   $L_{\infty}^s = 100$
#' 
#' - $K$ tells us how quickly we approach this asymptote as bears age.
#'  If we set $\mu = L_{\infty}/2$ and solve for $K$, this gives 
#'  $K = log(2)/(t-t_0)$. Then, noting that $\mu$ appears to be at 
#'  $L_{\infty}$/2 = 50 at around 3 years of age, we set 
#'  $K^s = log(2)/3 = 0.23$
#' 
#' -  $t_0$ is the x-intercept (i.e., the Age at which $\mu = 0$); 
#' an extrapolation of the curve to $x = 0$ suggests that $t_0$ 
#' should be slightly less than 0, say $t_0^s = -0.1$
#' 
#' - Most of the residuals when age is close to 0 are within 
#' approximately $\pm 10$ kg. For Normally distributed data, 
#' we expect 95% of the observations to be within $\pm 2\sigma$
#' so we set $\log(\sigma)^s = log(10/2) = 1.61$ 
#'
#' - We set $\theta^0$ = 1 (suggesting the variance increases 
#' linearly with the mean)
#'  
## ---------------------------------------------------------------
fitvb <- optim(par = c(100, 0.23, -0.1, 1.62, 1), 
               fn = logL,dat = beargrowth, 
               method = "BFGS", hessian = TRUE)
fitvb

#'  
#' Let's inspect the fit of our model by adding our predictions from the model to our data set (Figure \@ref(fig:vonb)).
#' 
## ----vonb, fig.cap="Fitted von Bertalanffy growth curve to weight-at-age data for black bears in Minnesota.", fig.align='center', out.width = "70%"----
pars<-data.frame(Linf = fitvb$par[1], K = fitvb$par[2], t0 = fitvb$par[3],
                 sigma = exp(fitvb$par[4]), theta = fitvb$par[5])
beargrowth <- beargrowth %>% 
  mutate(mu.hat = pars$Linf*(1 - exp(-pars$K*(age - pars$t0))))
ggplot(beargrowth, aes(age, wtkg)) + 
  geom_point() + geom_line(aes(age, mu.hat), col = "red", linewidth = 1) +
  ylab("Weight (kg)") + xlab("Age (years)") +
  theme_bw()


#' Next, let's create and inspect plots using standardized residuals 
## ----bresids, fig.cap="Standardized residuals (left panel) and square-rooted absolute values of these residuals (right panel) versus fitted values from the von Bertalanffy growth model applied to weight-at-age data from black bears in Minnesota.", out.width="90%", fig.align='center', fig.height=4, fig.width=8----
beargrowth <- beargrowth %>%
  mutate(sig.hats = pars$sigma^2*abs(mu.hat)^(2*pars$theta)) %>%
  mutate(stdresids = (wtkg - mu.hat)/(sig.hats),
         sqrt.abs.resids = sqrt(abs(stdresids)))
p1<-ggplot(beargrowth, aes(mu.hat, stdresids)) + geom_point() + 
  geom_hline(yintercept=0) + geom_smooth() + theme_bw()
p2<-ggplot(beargrowth, aes(mu.hat, sqrt.abs.resids)) + geom_point() + 
   geom_smooth() + theme_bw()
p1+p2

#' Variance covariance matrix of the parameters
(vcov.psi <- solve(fitvb$hessian))

#' 
#' The diagonals of this matrix tell us about the expected variability
#' of our parameter estimates (across repetitions of data collection 
#' and analysis), with the square root of these elements giving us 
#' our standard errors. Thus, we can form 95% confidence intervals 
#' for our model parameters using the assumption that the sampling 
#' distribution is Normal. For example, a 95% confidence interval for 
#' $L_{\infty}$ is given by:
#' 
## ---------------------------------------------------------------
SE.Linf<-sqrt(diag(vcov.psi))[1]
fitvb$par[1] + c(-1.96, 1.96)*rep(SE.Linf,2)

#' 
#' All good, but again, what if we want a confidence interval for 
#' $E[Y_i|Age_i] = \mu_i = L_{\infty}(1- e^{-K(Age_i-t_0)})$?
#' 
#' Options:
#' 
#' - Use a bootstrap (see Section \@ref(boot))
#' - Use the Delta method, which we will see shortly 
#' - Switch to Bayesian inference and use posterior distributions
#' 
#' ### Delta Method

## ---------------------------------------------------------------
age <- 1:35
mu.hat <- pars$Linf*(1-exp(-pars$K*(age-pars$t0)))
fprime<-as.matrix(cbind(1-exp(-pars$K*(age-pars$t0)),
                  pars$Linf*(age-pars$t0)*exp(-pars$K*(age-pars$t0)),
                  -pars$Linf*pars$K*exp(-pars$K*(age-pars$t0))))

#' 
#' We then determine the variance using matrix multiplication, pulling off the diagonal elements (the off-diagonal elements hold the covariances between observations for different ages):
## ---------------------------------------------------------------
var.mu.hat <- diag(fprime %*% solve(fitvb$hessian)[1:3,1:3] %*%t(fprime))

#' 
#' 
#' Alternatively, the `emdbook` package [@emdbook] has a function, `deltavar`, that will do the calculations for you if you supply a function for calculating $f()$ (via argument `meanval`) and you pass the asymptotic variance covariance matrix, $I^{-1}(\Psi)$ (via the `Sigma` argument).    
#' 
## ----warning=FALSE, message=FALSE-------------------------------
library(emdbook)
var.mu.hat.emd <- deltavar(linf*(1-exp(-k*(age-t0))), 
                      meanval=list(linf=fitvb$par[1], k=fitvb$par[2], t0=fitvb$par[3]), 
                      Sigma=solve(fitvb$hessian)[1:3,1:3])


#' We can take the square root of the diagonal elements of `var.mu.hat` 
#' to get SEs for forming pointwise confidence intervals for 
#' $E[Y_i|Age_i]$ at each Age.
#' 
## ----vbgci, fig.cap="Fitted von Bertalanffy growth curve and 95% confidence interval estimated from weight-at-age data for black bears in Minnesota."----
muhats <- data.frame(age=age, 
                     mu.hat=mu.hat, 
                     se.mu.hat = sqrt(var.mu.hat))
ggplot(beargrowth, aes(age, wtkg)) + 
  geom_point() + geom_line(aes(age, mu.hat), col="red", size=1) +
  ylab("Weight (kg)") + xlab("Age (years)") + theme_bw() +
  geom_ribbon(data=muhats, aes(x=age, 
                                 ymin=mu.hat-1.96*se.mu.hat, 
                                 ymax=mu.hat+1.96*se.mu.hat),
                 inherit.aes = FALSE, fill = "blue", alpha=0.5)

#' 
#' 
#' Lastly, we could calculate prediction intervals, similar to the example in Section \@ref(gls) by adding +/- 2$\sigma_i$, where $\sigma_i =\sqrt{\sigma^2 \mu_i^{2\theta}}$ (Figure \@ref(fig:piweigthage)). 
#' 
## ----piweigthage, fig.cap="Fitted von Bertalanffy growth curve and 95% confidence and prediction intervals estimated from weight-at-age data for black bears in Minnesota."----
muhats <- muhats %>%
  mutate(pi.up = mu.hat + 1.95*se.mu.hat + 2*sqrt(pars$sigma^2*mu.hat^(2*pars$theta)),
         pi.low = mu.hat - 1.95*se.mu.hat - 2*sqrt(pars$sigma^2*mu.hat^(2*pars$theta)))
ggplot(beargrowth, aes(age, wtkg)) + 
  geom_point() + geom_line(aes(age, mu.hat), col="red", size=1) +
  ylab("Weight (kg)") + xlab("Age (years)") + theme_bw() +
  geom_ribbon(data=muhats, aes(x=age, 
                                 ymin=mu.hat-1.96*se.mu.hat, 
                                 ymax=mu.hat+1.96*se.mu.hat),
                 inherit.aes = FALSE, fill = "blue", alpha=0.2)+
  geom_ribbon(data=muhats, aes(x=age, 
                                 ymin=pi.low, 
                                 ymax=pi.up),
                 inherit.aes = FALSE, fill = "red", alpha=0.2)

 