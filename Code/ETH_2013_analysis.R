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
cor(db1$N, db1$P) # maybe OK
# define inputs for translog
# define inputs for environmental variables
# define z inputs for regressing against technical
# efficiency

# core translog inputs
translog_inputs <- c("logN", "loglab")

# environmental inputs affecting production
# some variables are in log form because they
# have a lognormal distribution. Some are in
# a log form for computational ease
# we do not include legume because so few
# are present in the data
other_inputs <- c("noN", "logP", "logN*logP", "logslope",
                  "elevation", "irrig", "impr",
                   "rain_year", "SPEI",
                   "crop_count2", "logarea")

# variables explaining technical inefficiency.
# Note the sex is not included as virtually all
# households are headed by a male
# Also chose extension rather than ext_agent
# because overwhelmingly ext_agent == 1
# ed_any and literate likely measure the sa
zvars <- c("age", "ed_any",
           "micro_finance", "extension")

# basic form of the translog function
TL_form_basic <- translog_form("logyld", translog_inputs)

# 1. run a simple translog function and check
# for skewness and residuals histogram

TL_basic <- lm(TL_form_basic, data=db1)
hist(TL_basic$residuals, breaks=50) # skew left
skewness(residuals(TL_basic))

# 2. Incorporate the environmental variables
# and test again for skewness 

# TL
TL_form_full <- paste(TL_form_basic,
                      paste(other_inputs, collapse=" + "), sep=" + ")

TL_full <- lm(TL_form_full, data=db1)
hist(TL_full$residuals, breaks=50) # skew left but now less so
skewness(residuals(TL_full)) # skewed but now less so

# 3. skewness suggests the use of a SFA model

# translog estimation basic variables
sfa_TL_basic <- sfa(TL_form_basic,
                    data = db1)
sfa_res_basic <- summary(sfa_TL_basic)

# translog estimation full variables
sfa_TL_full <- sfa(TL_form_full,
                    data = db1)
sfa_res_full <- summary(sfa_TL_full)

# make a table of results
res_basic <- round(sfa_res_basic$mleParam, 3)
n <- row.names(res_basic)
p <- cut(res_basic[, 4], breaks=c(0, 0.01, 0.05, 1), labels=c("**", "*", ""), include.lowest = TRUE)
res_basic <- paste(res_basic[,1], " (", res_basic[, 2], ")", p, sep="")
res_basic <- data.frame(parameter=n, basic=res_basic)

res_full <- round(sfa_res_full$mleParam, 3)
n <- row.names(res_full)
p <- cut(res_full[, 4], breaks=c(0, 0.01, 0.05, 1), labels=c("**", "*", ""), include.lowest = TRUE)
res_full <- paste(res_full[,1], " (", res_full[, 2], ")", p, sep="")
res_full <- data.frame(parameter=n, full=res_full)

results <- full_join(res_basic, res_full)
move <- which(results$parameter %in% c("gamma", "sigmaSq") )
results <- rbind(results[-move, ], results[move,])

# table
knitr::kable(results)

# likelihood ratio test
lrtest(sfa_TL_full) # suggests sfa model fits better than OLS model

# likelihhod ratio suggests that full model fits better than basic



# 4. Use the z variables to find out what causes the gaps
# Look at barett paper which contains a good indictation
# of what are explanatory vars for this. The z vars are 
# regression the technical inefficiency term u againt
# explanatory factors, but still allowing for random noise
# in the estimation

TL_form_fullz <- paste(deparse(formula(TL_form_full), width.cutoff = 500),
                       paste(zvars, collapse=" + "), sep=" | ")

summary(sfa(TL_form_fullz, data=db1)) 

# 4 find maximums of functions -> look at Michiel's code for making a function to optimize
# once we have estimated the equation. We have all the components to numerically
# maximise the function

# use estimated frontier to get economic maximum
# compare actual of each farmer to the frontier
# Need a measure of potential from GYGA




# 5. endogeneity analysis
# as instruments we have several variables. We use dist_HH, age, years and other distance vars
# can also try the number of fert users in the local area
# and soil quality because of the relationship between the nitrogen and soil.



# note that we cannot include sex as it would appear 
# every single household in our sample is male
# literate and ed-any are highly correlated 80%
# probably still need to look for a better fit

library(AER)
m <- tobit(log(N) ~ log(dist_hh + 1) + log(dist_market + 1) +
                 irrig + impr + SPEI + manure + compost +legume +
                 log(slope+1) + elevation + 
                 extension + 
                 rain_year + age +
                 family_size + ed_any + relprice, data = db1) 

# battery of diagnostic plots for tobit model
yhat <- fitted(m)
rr <- resid(m, type = "response")
rd <- resid(m, type = "deviance")
apt <- m$y[, 1]

par(mfcol = c(2, 3))


plot(yhat, rr, main = "Fitted vs Residuals")
qqnorm(rr)
plot(yhat, rd, main = "Fitted vs deviance Residuals")
qqnorm(rd)
plot(apt, rd, main = "Actual vs deviance Residuals")
plot(apt, yhat, main = "Actual vs Fitted")

# now redo the sfa analysis, adding in the endogeneity part
# we know which values were left out of this stage of the analysis
indx <- m$na.action
db1.1 <- cbind(db1[-indx, ], rd)

# now run the translog again
TL_form_full_r <- paste(paste(deparse(formula(TL_form_full), width.cutoff = 500),
                      "rd", sep=" + "), "-noN", sep=" ")

sfa_TL_full_r <- sfa(TL_form_full_r,
                           data = db1.1)

# bootstrap the SEs -> problem with model matrix
# some terms do not appear
# bootstrapping should take about 25 minutes
# note that because we do not assume panel structure
# boostrapping within clusters is not required.
# But this would be required if we looked at
# a panel model.

B <- 100
coefM <- matrix(ncol=length(coef(sfa_TL_full_r)), nrow=B)
system.time({
for(i in 1:B){
  indx <- sample(1:nrow(db1.1), nrow(db1.1), replace=TRUE)
  modl <- sfa(formula(sfa_TL_full_r), data=db1.1[indx, ])
  coefM[i, ] <- coef(modl)
}
})
coef <- apply(coefM, 2, mean)
se <- apply(coefM, 2, sd)

N <- nrow(sfa_TL_full_r$dataTable)
df <- N - length(coef)
t <- coef/se
p.value = 2*pt(abs(t), df=df, lower=FALSE)

# Note that MPP cannot be calculated for plots with N=0 and are therefore set to 0.




# may need to take exponent to get back to
# highest point. Ask Michiel about tranforming values back
# and why they happen this way??
# coef <- coef(sfa_TL_full)
# logY <- fitted(sfa_TL_full)
# Y <- exp(logY)
# MPP <- (coef[2] + coef[5]*db1$logN + coef[18]*db1$logarea + coef[19]*db1$loglab + coef[21]*db1$logP) * Y/N
# MPPn <- ifelse(db1$N == 0, NA, 
#               exp(coef[2] + coef[5]*db1$logN + coef[18]*db1$logarea + coef[19]*db1$loglab + coef[21]*db1$logP) * (db1$yld/db1$N))


# work out economic yield using iterative algorithm as the
# place where the frontier meets the prices

# function for producing sum of frontier
# form to make calcualtions easier
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

# evaluate the production function on the frontier
modl <- sfa_TL_full
ysum <- Y_sum(modl)

# we only want the observations which are valid in
# the sfa estimation
indx <- modl$validObs
db1.1 <- db1[indx, ]

#x <- as.data.frame(modl$dataTable)
#names(x) <- gsub(":", "*", names(x))

logY <- with(db1.1, eval(parse(text=ysum)))
Y <- exp(logY)

# evaluate the MPP
MPP <- unlist(strsplit(ysum, "+", fixed=TRUE))
MPP <- MPP[grep("log(N)", MPP, fixed=TRUE)]
MPP <- gsub("log(N)*", "", MPP, fixed=TRUE)
MPP <- gsub("I(0.5 * log(N)^2)", "log(N)", MPP, fixed=TRUE)
MPP[1] <- gsub(" * log(N) ", "", MPP[1], fixed=TRUE)
# MPP <- gsub("logN", "log(N)", MPP)
MPP <- trimws(paste(MPP, collapse="+"))
MPP <- paste("(", MPP, ")", "*(Y/N) - 9.6", sep="")
# MPP <- ifelse(db1.1$N == 0, NA, with(db1.1, eval(parse(text=MPP))))

logarea <- quantile(db1.1$logarea)[4]
loglab <- quantile(db1.1$loglab)[4]
logP <- quantile(db1.1$logP)[4]

row=db1[5057, ]
row$N <- NULL

# function to optimize
MPP_f <- function(N){
  logY <- with(row, eval(parse(text=ysum)))
  Y <- exp(logY)
  # logarea <- row$logarea
  # loglab <- row$loglab
  # logP <- row$logP
  MPP = eval(parse(text=MPP))
  MPP
}



# need to find the optimum N using Newton-Raphson
# methods. But This would mean a whole bunch on N
# values


MPP_f <- function(N){

  
  # first get the frontier production function
  logY <- coef(TRA)[1]+coef(TRA)[2]*log(N) + coef(TRA)[3]*log(N)*log(N)
  Y = exp(logY)
  MPP = (coef(TRA)[2]+2*coef(TRA)[3]*log(N))*(Y/N)-avg_relprice # 9.6 is the N/maize price ratio.
  return(MPP)
}



# need to figure out whether to continue using noN in the equation
# What about a second stage tobit model? One that can be used to
# make a jusgement about the inefficiencies, especially if a number
# of them will be zero


# bootstrapping here because tobit is non-linear,
# but we do not need to worry about the households effects
# because we are using a cross section and have not
# assumed any structure in the data. This will
# give us the proper SEs



# residuals should be more negative if u term
# is more important

# 6. Analysing the actual yield gaps. Broken down
# into different terms. Need to get elasticities
# but these will be evaluated with other terms

