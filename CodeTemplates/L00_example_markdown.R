#'---
#'title: "Example"
#'author: ""
#'---
#' # A minimal example for knitr/ezknitr
#' 
#' This is an example of how to code (and then spin that code) to produce html documents with the ROxygen language
#' 
#' Regular, in-line text written like this with a "pound-appostrophe" symbol
#' 
#' ## Headers are coded with additional pounds
#' 
#' R-code is simply written without any pounds in "chunks," which have many options.
#' 
#' It's good practice to give each R-code chunk an informative name!
#' 
#' R-code comments are written as usual with one pound symbol
#' 
#' ## Document Preamble
#+ docPreamble, warning=FALSE, message=FALSE
library(knitr) 
# Chunk options
opts_chunk$set(fig.width=5, fig.height=5)

#' ## Data Entry
#+ dataEntry
# Snake dataset (Kery p 68)
mass <- c(6,8,5,7,9,11) # snake body mass in 10g
pop <- factor(c(1,1,2,2,3,3)) # population where factor = categorical variable
region <- factor(c(1,1,1,1,2,2)) 
hab <- factor(c(1,2,3,1,2,3)) # habitat
svl <- c(40, 45, 39, 50, 52, 57) # snout-vent length in mm?

#' ## Plotting
#' First, plot the regression between snout-vent length and body mass,
#' 
#' Then plot the residuals 
#+ massVsvl
# Plot of snout-vent length versus body mass
plot(svl,mass)
abline(lm(mass~svl))
# This plot will get saved as massVsvl-1.pdf with figure sizes specified above in the preamble
#+ massVsvlResiduals, fig.width=5, fig.height=4
# Plot of residuals from snout-vent length versus body mass
#+ fig.height=4, fig.width=8
plot(residuals(lm(mass~svl)), ylim = c(-1.5, 1.5))
abline(h = 0)
# This plot will get saved as massVsvlResiduals-1.pdf with figure sizes specified above in the chunk code

#' ## Document Footer
#'
#' Create a reproducible using:
#' 
#' rmarkdown::render("Rcode/example_roxygen.R", output_dir = "Output", output_format="pdf_document")
#' 
#' Alternative output formats include "html_document", "word_document",or "all" 
#'
#' - Also, notice how easily you can add bullet points with roxygen!
#' 
#' Session Information:
#+ sessionInfo
sessionInfo()