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
library(moments)
library(dplyr)

# source in prepared data and functions
source(file.path(root, "Code/ETH_2013_prepare4analysis.R"))
source(file.path(root, "Code/translog_formula.R"))
source(file.path(root, "Code/sfaFormEval.R"))

#' ------------------------------------------------------------------------------------------------
#' Define inputs for use in translog functions
#' ------------------------------------------------------------------------------------------------

# core translog inputs
translog_inputs <- c("logN", "loglab")

# environmental inputs affecting production
other_inputs <- c("noN", "logslope",
                  "elevation", "irrig", "impr",
                  "rain_year", "SPEI",
                  "crop_count2", "logarea")

# variables explaining technical inefficiency.
zvars <- c("age", "ed_any",
           "micro_finance", "extension")

#' ------------------------------------------------------------------------------------------------
#' Define translog formulas
#' ------------------------------------------------------------------------------------------------

# basic specification
TL_form_basic <- translog_form("logyld", translog_inputs)

# full specification including environmental inputs
TL_form_full <- paste(TL_form_basic,
                      paste(other_inputs, collapse=" + "), sep=" + ")

# full specification including z variables
TL_form_fullz <- paste(deparse(formula(TL_form_full), width.cutoff = 500),
                       paste(zvars, collapse=" + "), sep=" | ")

#' ------------------------------------------------------------------------------------------------
#' Estimate OLS models and check for skewness
#' Of cobb douglass and translog functions
#' ------------------------------------------------------------------------------------------------

# Cobb Douglass



# Translog



#' ------------------------------------------------------------------------------------------------
#' Run Tobit model, first stage endogeneity analysis
#' ------------------------------------------------------------------------------------------------




#' ------------------------------------------------------------------------------------------------
#' Incorporate tobit residuals into analysis for
#' second stage endogeneity analysis with bootstrapping
#' ------------------------------------------------------------------------------------------------



#' ------------------------------------------------------------------------------------------------
#' Find MPP and optimal Nitrogen. This cannot be done
#' analytically for a translog => numerical methods
#' ------------------------------------------------------------------------------------------------

# model that we want to find optimum for
modl <- sfa(TL_form_full, data=db1)

# function to evalute optimum
Nopt_f <- function(N, modl, row){
  ysum <- sfaFormEval(modl)
  ysum <- gsub("logN", "log(N)", ysum)
  logY <- with(row, eval(parse(text=ysum)))
  Y <- exp(logY)
  MPP = with(row, ((coef(modl)["logN"] + 
            2*coef(modl)["I(logN^2)"]*log(N) +
            coef(modl)["logN:loglab"]*loglab)*
           (Y/N)) -  relprice)
  return(MPP)
}



# we can only look at those plots which have
# been used in estimating modl. In addition, we only
# want those plots that have N > 0

db1.1 <- db1[modl$validObs, ] # removes NA, and NaN values
db1.1 <- db1.1[db1.1$N > 0, ] # removes rows for which N was zero
Nobs <- db1.1$N # save N values for later
db1.1$N <- NULL


# function to iterate over. Finds optimal N 
# for each plot in the data. However, in order
# to use the uniroot function for this, we need
# to make sure that the lower and upper
# limit of the interval argument evaluate to
# values with opposite signs. Hence the while 
# and if statements in the following function
f <- function(i){
  row <- db1.1[i, ]
  x <- 10
  while(x > 1){
    if(Nopt_f(x, modl, row) < 0){
      x <- x - 1
    }else{
      break}
  }
  
  if(x == 1){
    return(NA)
    }else{
  root <- uniroot(function(x) {Nopt_f(x, modl, row)}, interval=c(x, 2000))$root
  return(root)
  }
}

# iterate over each plot
db1.1$Npm <- sapply(1:nrow(db1.1), f)
db1.1$N <- Nobs
db1.1$Ndif = db1.1$N - db1.1$Npm


#' ------------------------------------------------------------------------------------------------
#' MPP
#' ------------------------------------------------------------------------------------------------

MPP_f <- function(row){
  ysum <- sfaFormEval(modl)
  ysum <- gsub("logN", "log(N)", ysum)
  logY <- with(row, eval(parse(text=ysum)))
  Y <- exp(logY)
  MPP = with(row, ((coef(modl)["logN"] + 
                      2*coef(modl)["I(logN^2)"]*log(N) +
                      coef(modl)["logN:loglab"]*loglab)*
                     (Y/N)))
  return(MPP)
}

f2 <- function(i){
  row <- db1.1[i, ]
  MPP_f(row)
}

db1.1$MPP <- sapply(1:nrow(db1.1), f2)

#' ------------------------------------------------------------------------------------------------
#' calculate the yield gaps
#' ------------------------------------------------------------------------------------------------

sumzone_sfaTL <- db1.1 %>% group_by(ZONE) %>%
  summarize(
    N=mean(N, na.rm=T),
    Npm=mean(Npm, na.rm=T),
    MPPmean=mean(MPP[!is.infinite(MPP)], na.rm=T),
    MVC=mean((Pm[!is.infinite(MPP)]*MPP[!is.infinite(MPP)])/Pn[!is.infinite(MPP)], na.rm=T),
    Ndif=mean(Ndif, na.rm=T),
    Number=n())