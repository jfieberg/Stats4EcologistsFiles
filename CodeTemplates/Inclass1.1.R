#' ---
#' title: "L01: Lion Nose example from class"
#' author: "(leave blank for HW assingments)"
#' date: ""
#'---

#' # InClass1.1: Simulations to understand sampling distributions
#' 
#' Includes:
#' 
#' 1. Lion noses linear regression
#' 2. Data generation consistent with model
#' 3. Linear regression of this first dataset
#' 4. In-class Sampling Distribution Simulation Assignment
#' 
#' ## Document Preamble
#' 
#' Load libraries
#+ Warning=FALSE, message=FALSE
library(knitr)
library(abd)

#' Settings for Knitr (optional)
#+ simsKnitSettings
opts_chunk$set(fig.width = 8, fig.height = 6)
set.seed(0204)

#' ## 1. Lion noses linear regression: 
#' 
#' ####  Data entry
data(LionNoses)
head(LionNoses)

#' #### Fit linear model
lm.nose<-lm(age~proportion.black, data=LionNoses)

#' #### Parameters:
 
#' Coefficients and residual variation are stored in lmfit:
coef(lm.nose)
summary(lm.nose)$sigma # residual variation

#' What else is stored in lmfit? (residuals, variance covariance matrix, etc)
names(lm.nose)
names(summary(lm.nose))

#' ## 2. Data generation consistent with fitted model

## Use the same sampmle size Sample size - use length so it matches sample size of original data
n <- length(LionNoses$age)

## Predictor - copy of original proporation black data, now in vector
p.black <- LionNoses$proportion.black

## Parameters
#sigma <- summary(lm.nose)$sigma # residual variation
#betas <- coef(lm.nose)# regression coefficients
betas<-c(0.88, 10.65)# regression coefficients
sigma<- 1.67 # residual variation

## Errors and response
# Residual errors are modeled as ~ N(0, sigma)
epsilon <- rnorm(n, 0, sigma)

# Response is modeled as linear function plus residual errors
y <- betas[1] + betas[2]*p.black + epsilon

#' ## 3. Linear regression of this generated dataset

# Fit of model to simulated data:  
lmfit.generated <- lm(y ~ p.black)
summary(lmfit.generated)

#' ----------------------------------------------------------------------------

#' ## In-Class Sampling Distribution Simulation Assignment 
#'

#' 
#' ### Exercise 1: 
#' 
#' 1. Generate 5000 datasets using the same code
#' 2. Fit a linear regression model to each dataset "lm.temp"
#' 3. Store the estimates of $\beta_1$
#' 
#' Hint: if you get stuck, try starting with a small number of simulations (less than 5000) until you get the code right. 

#	set up a matrix of size 1 by 5000 to store our estimates of beta_1
beta.hat<- matrix(NA,	nrow	=	5000,	ncol	=	1)

# For loop: 
   
for(i	in	1:5000){
  # Insert	R	code here to simulate data, fit the model, and store 
  # results
}

#' ----------------------------------------------------------------------------
#' 
#' ### Exercise 2: 

#' Repeat the simulation exercise, but this time, also calculate a
#' t-statistic for each model fit
#' 
#	set up a matrices to store results here  

# For loop: 

for(i	in	1:5000){
  # Insert	R	code here to simulate data, fit the model, calculate 
  # t statistic, and store results
  
}




#' ----------------------------------------------------------------------------
#' 
#' ### Exercise 3: 

#' Now, add a line of code that calculates a confidence interval for
#' $\beta_1$ each time you fit a model to a newly simulated data set.
#' 
#	set up a matrices to store results here  

# For loop: 

for(i	in	1:5000){
  # Insert	R	code here to simulate data, fit the model, 
  # calculate t statistic, and store results
  
}


#'   
#' ## Document footer
#' 
#' Reproducible report created using:
#' 
#' rmarkdown::render("  ", output_dir = "  ", output_format="pdf_document")
#' 
#' Alternative output formats include "html_document", "word_document",or "all" 
#' 


#'
#' ### Session Information:
sessionInfo()