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
library(qdap)
library(ggplot2)
library(AER)

# source in prepared data and functions
# source(file.path(root, "Code/ETH_2013_prepare4analysis.R"))
db1 <- readRDS("Cache/db1.rds")
source(file.path(root, "Code/translog_formula.R"))
source(file.path(root, "Code/sfaFormEval.R"))

#' ------------------------------------------------------------------------------------------------
#' Define inputs for use in translog functions
#' ------------------------------------------------------------------------------------------------

# core translog inputs
translog_inputs <- c("log(N)", "log(lab)")

# environmental inputs affecting production
environmental <- c("log(slope)",
                  "elevation",
                  "log(area)",
                  "ph", "ph2", "SOC")

# environmental and all dummy variables
dummies <- c("noN", "irrig", "impr", "crop_count2")

# variables explaining technical inefficiency.
zvars <- c("-1", "age", "ed_any",
           "micro_finance", "extension",
           "log(off_farm_income)", "log(area_tot)",
           "log(HHEA)", "log(dist_market)")

# for some variables that we take the log of, we need
# to add 1
db1$N <- db1$N + 1
db1$lab <- db1$lab + 1
db1$P <- db1$P + 1
db1$slope <- db1$slope + 1
db1$off_farm_income <- db1$off_farm_income + 1



#' ------------------------------------------------------------------------------------------------
#' Define translog formulas
#' ------------------------------------------------------------------------------------------------

# basic specification
TL_form_basic <- translog_form("log(yld)", translog_inputs)

# full specification including environmental inputs
TL_form_full <- paste(TL_form_basic,
                      paste(environmental, collapse=" + "),
                      paste(dummies, collapse=" + "), sep=" + ")

# full specification including z variables
TL_form_fullz <- paste(deparse(formula(TL_form_full), width.cutoff = 500),
                       paste(zvars, collapse=" + "), sep=" | ")

#' ------------------------------------------------------------------------------------------------
#' Estimate OLS models and check for skewness
#' Of cobb douglass and translog functions
#' ------------------------------------------------------------------------------------------------

# Cobb Douglass
CD1 <- lm(log(yld) ~ log(N) + log(lab), data=db1)
CD2_form <- paste(deparse(formula(CD1), width.cutoff = 500),
                  paste(environmental, collapse=" + "),
                  paste(dummies, collapse=" + "), 
                  sep=" + ")
CD2 <- lm(CD2_form, data=db1)

par(mfrow=c(1, 2))
hist(residuals(CD1), main="residuals CD core model")
hist(residuals(CD2), main="residuals CD full model")

# Translog
TL1 <- lm(TL_form_basic, data=db1)
TL2 <- lm(TL_form_full, data=db1)

par(mfrow=c(1, 2))
hist(residuals(TL1), main="residuals CD core model")
hist(residuals(TL2), main="residuals CD full model")

#' ------------------------------------------------------------------------------------------------
#' Run sfa models
#' ------------------------------------------------------------------------------------------------

# Cobb Douglass
sfaCD1 <- sfa(log(yld) ~ log(N) + log(lab), data=db1)
sfaCD2 <- sfa(CD2_form, data=db1)

par(mfrow=c(1, 2))
hist(residuals(CD1), main="residuals CD core model")
hist(residuals(CD2), main="residuals CD full model")

# Translog
sfaTL1 <- sfa(TL_form_basic, data=db1)
sfaTL2 <- sfa(TL_form_full, data=db1)

# lrtest to determine what is a better fit


#' ------------------------------------------------------------------------------------------------
#' z variable analysis
#' ------------------------------------------------------------------------------------------------

sfaTL_z <- sfa(TL_form_fullz, data=db1)


#' ------------------------------------------------------------------------------------------------
#' Run Tobit model, first stage endogeneity analysis
#' ------------------------------------------------------------------------------------------------

m <- tobit(log(N) ~ relprice + log(dist_market), data = db1)

log(dist_market) + family_size
             ed_any + age + sex +
             extension  + impr + 
             log(slope+1) + elevation + 
             AI + log(TS) + ph + SOC
              twi,

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

#' ------------------------------------------------------------------------------------------------
#' Incorporate tobit residuals into analysis for
#' second stage endogeneity analysis with bootstrapping
#' ------------------------------------------------------------------------------------------------

indx <- m$na.action
db1.1 <- cbind(db1[-indx, ], rd) # use the same 

# now run the translog again
TL_form_full_r <- paste(paste(deparse(formula(TL_form_full), width.cutoff = 500),
                              "rd", sep=" + "), "-noN", sep=" ")

sfa_TL_full_r <- sfa(TL_form_full_r,
                     data = db1.1)

# standard errors with a nonlinear (tobit) function
# are no longer valid -> as a result we have to 
# bootstrap the SEs to get valid p-values

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


#' ------------------------------------------------------------------------------------------------
#' Find MPP and optimal Nitrogen. This cannot be done
#' analytically for a translog => numerical methods
#' ------------------------------------------------------------------------------------------------

# if there is no relative price we cannot
# work out the optimum
db1 <- db1[!is.na(db1$relprice),] # remove NA from relprice

# model that we want to find optimum for
modl <- sfa(TL_form_full, data=db1)

# function to evalute optimum
Nopt_f <- function(N, modl, row){
  ysum <- sfaFormEval(modl)
  logY <- with(row, eval(parse(text=ysum)))
  Y <- exp(logY)
  MPP = with(row, ((coef(modl)["log(N)"] + 
            2*coef(modl)["I(log(N)^2)"]*log(N) +
            coef(modl)["log(N):log(lab)"]*log(lab))*
           (Y/N)) -  relprice)
  return(MPP)
}


# we can only look at those plots which have
# been used in estimating modl. In addition, we only
# want those plots that have N > 0
db1.1 <- db1[modl$validObs, ] # removes NA, and NaN values that cannot be used in estimation
Nobs <- db1.1$N # save N observed values for later
db1.1$N <- NULL # remove N values until after optimum is found


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
  print(i)
  while(x > 1){
    if(Nopt_f(x, modl, row) < 0){
      x <- x - 1
    }else{
      break}
  }
  
  if(x == 1){
    return(NA)
    }else{
  root <- uniroot(function(x) {Nopt_f(x, modl, row)}, interval=c(x, 4000))$root
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

# Function for calcualting the MPP
MPP_f <- function(row){
  ysum <- sfaFormEval(modl)
  logY <- with(row, eval(parse(text=ysum)))
  Y <- exp(logY)
  MPP = with(row, ((coef(modl)["log(N)"] + 
                      2*coef(modl)["I(log(N)^2)"]*log(N) +
                      coef(modl)["log(N):log(lab)"]*log(lab))*
                     (Y/N)))
  return(MPP)
}

# function to apply MPP to each plot
f2 <- function(i){
  row <- db1.1[i, ]
  if(row$N == 1){
    return(NA)
  } else {
  return(MPP_f(row))
  }
}

# get MPP values
db1.1$MPP <- sapply(1:nrow(db1.1), f2)

#' ------------------------------------------------------------------------------------------------
#' calculate the four types of yield
#' 1. Technical Efficiency Yield
#' 2. Economic Yield
#' 3. Feasible Yield
#' 4. Potential Yield
#' ------------------------------------------------------------------------------------------------

# frontier translog model
model <- sfa(TL_form_full, data=db1.1)

# 1. Technical efficiency yield is found using the
# output of the sfa model
db1.1 <- db1.1 %>%
  rename(Y = yld) %>%
  mutate(
    Ycor = exp(as.numeric(fitted(model))+log(as.numeric(efficiencies(model)))), 
    err = Ycor-Y,
    TEY = exp(as.numeric(fitted(model))),
    TE = as.numeric(efficiencies(model)),
    resid = as.numeric(resid(model))
  )

table((db1.1$Y - db1.1$TEY) > 0) # 334

# 2. Economic yield is found by evaluating the frontier 
# function at the economically optimal nitrogen rate (Npm)

# introduce a cap on nitrogen
# (check this from experimental plots of Ethiopia)
Npy <- 300

# Cap Npm
db1.1 <- mutate(db1.1, Npm = ifelse(Npm > Npy, Npy, Npm))

# economic yield is evaluated at the frontier with
# optimal nitrogen but everything else stays the same

ysum <- sfaFormEval(modl)
ysum <- mgsub("log(N)", "log(Npm)", ysum)
logY <- with(db1.1, eval(parse(text=ysum)))
db1.1$EY <- exp(logY)

# 3. PFY: Feasible yield
# To improve this part, we could also argue that: (1) hybrid seeds are used, (2) pestices are used, (3) higher levels of capital and labour are used.
# We assume that all farmers use pesticides and increase assets and labour with 10%

# CHECK increase seed QUANTITY
# what should be multiplied by 1.1 here? loglab of log(lab*1.1)??
library(qdap)
ysum <- sfaFormEval(modl)
ysum <- mgsub(dummies,"1",ysum)
ysum <- mgsub(c("lab"), paste0(c("lab"), "*1.1"), ysum)
ysum <- mgsub("log(N)", "log(Npy)", ysum)

db1.1$Npy <- Npy
db1.1$PFY <-  exp(with(db1.1, eval(parse(text=ysum))))


dbgap <- select(db1.1, hhid, ea_id, parcel_id, field_id, AEZ,
                Y, Ycor, TEY, EY, PFY)


# 4. PY: Potential yield and HFY: Highest farmers yield
# Merge Yield potential with maize plot database
db1.2 <- db1 %>% select(hhid, parcel_id, field_id, PY = YW) %>% 
  na.omit %>% unique() %>%
  mutate(PY=PY*1000) %>% 
  left_join(dbgap,.) %>%
  do(filter(., complete.cases(.))) %>%
  ungroup() %>%
  group_by(AEZ) %>%
  mutate(HFY90 = quantile(Y, probs = 0.9, na.rm = TRUE),
         HFY95 = quantile(Y, probs = 0.95, na.rm = TRUE),
         HFY100 = quantile(Y, probs = 1, na.rm = TRUE)) %>%
  ungroup()

#' ------------------------------------------------------------------------------------------------
#' Yield gap consistency checks and corrections - speak to Michiel
#' ------------------------------------------------------------------------------------------------

#  We cap all values at PY because we consider this as an absolute potential and recalculate all gaps.
db1.3 <- mutate(db1.2, PFY = ifelse(PY-PFY<0, PY, PFY),
              EY = ifelse(PY-EY<0, PY, EY),
              TEY = ifelse(PY-TEY<0, PY, TEY),
              Ycor = ifelse(PY-Ycor<0, PY, Ycor),
              Y = ifelse(PY-Y<0, PY, Y))

# Calculate TYG using UY as reference
db1.4 <- db1.3 %>% 
  mutate(
    ERROR_l = Ycor-Y,      # Error gap
    ERROR_s = Y/Ycor,      # Error gap
    TEYG_l = TEY-Ycor,     # Technical efficiency yield gap using Ycor as basis
    TEYG_s = Ycor/TEY,     # Technical efficiency yield gap using Ycor as basis
    EYG_l = EY-TEY,        # Economic yield gap
    EYG_s = TEY/EY,        # Economic yield gap
    EUYG_l = PFY-EY,       # Feasible yield gap
    EUYG_s = EY/PFY,       # Feasible yield gap
    TYG_l = PY-PFY,        # Technology yield gap
    TYG_s = PFY/PY,        # Technology yield gap
    YG_l = PY-Y,           # Yield gap
    YG_s = Y/PY,           # Yield gap
    YG_l_Ycor = PY-Ycor,   # Yield gap with Ycor as reference
    YG_s_Ycor = Ycor/PY)   # Yield gap with Ycor as reference


# Check if separate yield gaps add up to total yield gap
X_overall <- db1.4 %>%
  mutate(check_l = YG_l/(ERROR_l + TEYG_l + EYG_l + EUYG_l + TYG_l), # Note that for a small number of observatios YG_l=0 resulting in 0/0 which is NaN
         check_s = YG_s/(ERROR_s * TEYG_s * EYG_s * EUYG_s * TYG_s),
         check_l2 = YG_l_Ycor/(TEYG_l + EYG_l + EUYG_l + TYG_l),
         check_s2 = YG_s_Ycor/(TEYG_s * EYG_s * EUYG_s * TYG_s))
summary(X_overall)


# Create database with relevant variables for further analysis
db1.5 <- select(db1.4, hhid, ea_id, parcel_id, field_id, AEZ,Y,
                          Ycor, TEY, EY, PFY, HFY90, HFY95, HFY100, PY, ERROR_l, ERROR_s, TEYG_l, TEYG_s, EYG_l, EYG_s, 
                          EUYG_l, EUYG_s, TYG_l, TYG_s, YG_l, YG_s, YG_l_Ycor, YG_s_Ycor)

# Clean up
rm(list =grep("X_", ls(), value = TRUE)) # remove checks
rm(db1, db1.1, db1.2, db1.3, db1.4)

# Save data
saveRDS(db1.5, "Cache/db1.5.rds")


# sumzone_sfaTL <- db1.1 %>% group_by(ZONE) %>%
#   summarize(
#     N=mean(N, na.rm=T),
#     Npm=mean(Npm, na.rm=T),
#     MPPmean=mean(MPP[!is.infinite(MPP)], na.rm=T),
#     MVC=mean((Pm[!is.infinite(MPP)]*MPP[!is.infinite(MPP)])/Pn[!is.infinite(MPP)], na.rm=T),
#     Ndif=mean(Ndif, na.rm=T),
#     Number=n())
