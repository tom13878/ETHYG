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

# define inputs for translog
# define inputs for environmental variables
# define z inputs for regressing against technical
# efficiency

# core translog inputs
translog_inputs <- c("logN", "logarea", "loglab")


other_inputs <- c("noN", "logP", "logN*logP", "logslope",
                  "elevation", "irrig", "impr",
                  "rain_year", "SPEI",
                  "crop_count2")

TL_form_basic <- translog_form("logyld", translog_inputs)
TL_form_full <- paste(TL_form_basic,
                      paste(other_inputs, collapse=" + "), sep=" + ")

# run sfa model
modl <- sfa(TL_form_basic, data=db1)

# helper function for creating sum
Y_sum <- function(modl){
  
  # get the form of the frontier equation
  form <- strsplit(as.character(formula(modl)), "~ ")[[1]][2]
  
  # combine coefs with variables
  form <- unlist(strsplit(form, " + ", fixed=TRUE))
  coef <- coef(modl)[2:(length(coef(modl))-2)]
  names(coef) <- gsub(":", "*", names(coef))
  coef <- coef[form]
  
  # paste together
  form <- paste(paste(coef,
                      form, sep = " * "), collapse=" + ")
  paste(coef(modl)[1], form, sep=" + ")
}

ysum <- Y_sum(modl)
logY <- with(db1, eval(parse(text=ysum)))
Y <- exp(logY)
ysum <- gsub("logN", "log(N)", ysum)


row <- db1[4412, ]
row$N <- NULL

# MPP
MPP_f <- function(N){
  logY <- with(row, eval(parse(text=ysum)))
  Y <- exp(logY)
  logarea <- row$logarea
  loglab <- row$loglab
  logP <- row$logP
  MPP = (coef(modl)[2] + coef(modl)[5] * log(N) + coef(modl)[8] * logarea + coef(modl)[9] * loglab) * (Y/N)
  MPP
}

uniroot(MPP_f, interval=c(10, 1000))

N <- seq(1, 100, length.out=100)
plot(N, MPP_f(N))


