#'=================================================================================================
#' Project:  IMAGINE ETH
#' Subject:  Analysis file
#' Author:   Michiel van Dijk & Tom Morley
#' Contact:  michiel.vandijk@wur.nl, Tomas.morley@wur.nl
#' Output:   Results
#'=================================================================================================

# set project root
library(rprojroot)
root <- find_root(is_rstudio_project)

# load packages
library(frontier)
library(stargazer)

# source in prepared data and functions
source(file.path(root, "Code/ETH_2013_prepare4analysis.R"))
source(file.path(root, "Code/translog_formula.R"))
source(file.path(root, "Code/sfaTable.R"))

# simple stochastic frontiers analysis with a cobb douglass
# production function introduce types of variables one at a
# time. Look in henningsen for ways to compare the translog and
# the cobb  douglass functions. Because the cobb douglass is
# nested in the translog we can use a likelihood ratio test
# to ocmpare the two types of model
# lrtest: likelihood ratio test function
# then look at skewness of the residuals and plot histograms

# 1. Core production variables and those of scientific validity

# core variables
inputs <- c("logN", "logarea", "loglab")

# cobb douglass estimation
sfaCD1 <- sfa(logyld ~ logN + loglab + logarea,
              data = db1)

# make a plot of the frontier. Every variable
# except for Nitrogen is evaluated at its 75th
# percentile
X <- model.matrix(logyld ~ logN + loglab + logarea, data=db1)
Xnew <- data.frame(`/(Intercept/)` = 1,
                   `logN` = X[,2],
                   `loglab` = quantile(X[, 3])[4],
                   `logarea` = quantile(X[, 4])[4])
ypred <- Xnew %*% coef(sfaCD1)[1:ncol(X)]
plot(exp(X[,2]), exp(ypred))

# form of the translog
trans1 <- translog_form("logyld", c("logN", "loglab", "logarea"))

a = 10
eval(parse(text="a+5"))
ypred <- with(Xnew,
eval(parse(text=strsplit(trans1, "~ ")[[1]][2])))

# combine coefs with variables
test <- strsplit(trans1, "~ ")[[1]][2]
test2 <- unlist(strsplit(test, " + ", fixed=TRUE))
test3 <- paste(paste(coef(sfaTL1)[-c(1,11, 12)], test2, sep = " * "), collapse=" + ")

# calcualte predicted value, still need to add the intercept. but on the right track
ypred <- with(Xnew,
              eval(parse(text=test3)))

plot(exp(X[,2]), exp(ypred))

# translog function
sfaTL1 <- sfa(trans1,
     data = db1)
# use the formula to carry out an actual sum

ypred <- Xnew %*% coef(sfaCD1)[1:ncol(X)]
plot(exp(X[,2]), exp(ypred))

waldtest(sfaCD1, sfaTL1) # reject cobb douglass in favour of translog
lrtest(sfaCD1, sfaTL1) # reject cobb douglas in favour of translog

# 2. expand the selection to include the environmental
# variables as well, these include the plot level
# characteristics and all the climatic variables

sfaCD1 <- sfa(logyld ~ logN + loglab + logarea,
              data = db1)

# 3.  

# 4. endogeneity 


# 4. 

