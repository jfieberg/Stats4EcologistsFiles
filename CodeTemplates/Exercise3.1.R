#' ---
#' title: "Homework 1. Linear Regression"
#' author: "Leave Blank for grading"
#' date: ""
#'---
#' 

#' ## Document Preamble
#+ warning=FALSE, message=FALSE
# Load Libraries
 library(knitr)
 library(here)

# Set knitr options
opts_chunk$set(fig.width = 6, fig.height=5)

# Clear Environment (optional)
remove(list=ls())

# Set seed 
set.seed(314159)

#' ## Section 1. Bootstrapping RIKZ data
#' 
#' Data Entry
# Read in data from .csv and look at first 6 rows
rikzData <- read.csv(here("data", "RIKZdat.csv"))
head(rikzData)

#' Linear regression analysis
#+ linReg01
# Simple linear regression and summary
lmfit <- lm(Richness~exposure, data=rikzData)
summary(lmfit)

# Simple linear regression confidence intervals for intercept/slope
confint(lmfit)

#' Cluster-level bootstrap example in class (one bootstrap)
# Data processing
uid <- unique(rikzData$Beach) # unique id for each beach
nBeach <- length(uid) # number of beaches

### One bootstrap:
# Take a sample from x (uid) of size nBeach with replacement:
bootIDs <- data.frame(Beach = sample(x = uid, size = nBeach, replace = TRUE))
bootIDs

# Use this to sample from original data by beach
bootDat <- merge(bootIDs, rikzData)
table(bootDat$Beach) ## this table shows how many obs are drawn for each beach in bootstrap sample.

# Double check sample sizes worked (these should match):
length(rikzData$Beach) # original data
length(bootDat$Beach) # bootstrap sample

#' Complete cluster-level bootstrap for nonindependence

# Set up object to store results
# 


# Loop


# Plot distribution of bootstrap estimates
 

# Bootstrap confidence intervals (e.g., using quantile function)
 

# Compare to linear regression confidence intervals:
confint(lmfit)

#' ### Interpretation
#'
  

#' ## Section 2. Kery Ch 6.3 and 6.4
#' 

#' ## Document footer 
#' 
#' NOTE: feel free to change the name of the file containing your solution (and modify the line of code, below)
#' 
#' Reproducible report created using:
#' 
#' rmarkdown::render("homework/HW01_Template.R", output_dir = "homework/output", output_format="pdf_document")
#' 
#' ### Session Information:
#+ sessionInfo
sessionInfo()