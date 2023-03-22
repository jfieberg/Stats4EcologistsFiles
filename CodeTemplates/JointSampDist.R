#' ## Multivariate Sampling Distribution
#' 
#'
#' True parameterse
betas <- c(0.5, 1) 
sigma<-1.5

#' Simulation setup
nsims <- 5000 # number of simulations
n<-30
beta.hat<- matrix(NA,   nrow    =   nsims,  ncol    =   2)
x<-runif(30, 0, 1) # single predictor

# Simulation
for(i in 1:nsims){
  epsilon <- rnorm(n, 0, sigma) # random errors
  y <- betas[1] + betas[2]*x + epsilon # response
  lm.temp <- lm(y ~ x)
  ## extract beta-hat  
  beta.hat[i,] <- coef(lm.temp) 
}

#' The slope and the intercept walk in the room together, and when we tend to
#' get a large intercept, we tend to get a smaller slope.  The estimates of
#' the intercept and slope are negatively correlated.
cor(beta.hat)
plot(beta.hat[,1], beta.hat[,2], 
     xlab=expression(hat(beta[0])), 
     ylab=expression(hat(beta[1])))
Sampdist <- ks::kde(x = beta.hat) 
plot(Sampdist, display = "persp",
  col.fun = viridis::viridis,
  xlab = expression(hat(beta[0])),
  ylab = expression(hat(beta[1]))
) 


#' So...we need to think of their **joint sampling distribution**.


#' When our sample size is large, the sampling distribution of MLEs is given by
#' $(\hat{\beta}) \sim N(\beta, \Sigma)$ where $\Sigma$ is a variance/covariance 
#' matrix with the sqrt of the diagonals giving us the SEs for each parameter
#' and off diagonals describing how parameters co-vary across different samples. 

#' The diagonals of $\Sigma$ give us the variances of $\hat{\beta}$ and the 
#' off diagonals tell us how $\hat{\beta_1}$ and $\hat{\beta}_1$ covary.

