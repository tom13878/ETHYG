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
library(stargazer)

# source in prepared data and functions
# source(file.path(root, "Code/ETH_2013_prepare4analysis.R"))
db1 <- readRDS("Cache/db1.rds")
db1 <- unique(db1)
source(file.path(root, "Code/translog_formula.R"))
source(file.path(root, "Code/sfaFormEval.R"))
source(file.path(root, "Code/sfa_tab.R"))

# labour values this high are ridiculous
# CHECK NEED TO WINSOR ALL VALUES I THINK
db1 <- db1[db1$lab < 3000, ]

#' ------------------------------------------------------------------------------------------------
#' Define inputs for use in translog functions
#' ------------------------------------------------------------------------------------------------

# core translog inputs
translog_inputs <- c("log(N)", "log(lab)")

# environmental inputs affecting production
environmental <- c("log(slope)",
                  "elevation",
                  "log(area)",
                  "SOC",
                  "log(rain_wq)")

# GYGA variables
GYGA <- c("GGD", "AI", "TS")

# environmental and all dummy variables
dummies <- c("noN", "impr", "crop_count2", "phdum_gt70", "phdum55_2_70")

# variables explaining technical inefficiency.
zvars <- c("-1", "extension", "age",
           "I(age^2)", "ed_any", "sex", "title",
           "log(area_tot)")


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
TL_form_environ <- paste(TL_form_basic,
                      paste(environmental, collapse=" + "),
                      paste(dummies, collapse=" + "), sep=" + ")

# spec with GYGA variables
TL_form_environ_GYGA <- paste(TL_form_environ,
                              paste(GYGA, collapse=" + "),
                               sep=" + ")

# full specification including z variables
TL_form_fullz <- paste(deparse(formula(TL_form_environ_GYGA), width.cutoff = 500),
                       paste(zvars, collapse=" + "), sep=" | ")

TL_form_fullzr <- paste(paste(deparse(formula(TL_form_environ_GYGA), width.cutoff = 500), "rd - noN", sep=" + "),
                       paste(zvars, collapse=" + "), sep=" | ")

#' ------------------------------------------------------------------------------------------------
#' Summary statistics
#' ------------------------------------------------------------------------------------------------

# summary stats of some of the important variables
vars <- c("N", "lab", "area", "slope", "elevation",
          "SOC", "rain_wq", "GGD", "AI", "TS",
          "yesN", "impr", "extension", "title")
dbsum <- db1[, vars]

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

# par(mfrow=c(1, 2))
# hist(residuals(CD1), main="residuals CD core model")
# hist(residuals(CD2), main="residuals CD full model")

# Translog
TL1 <- lm(TL_form_basic, data=db1)
TL2 <- lm(TL_form_environ, data=db1)
TL3 <- lm(TL_form_environ_GYGA, data=db1)

stargazer(CD1, CD2, TL1, TL2, TL3, type = "text")

# par(mfrow=c(1, 2))
# hist(residuals(TL1), main="residuals CD core model")
# hist(residuals(TL2), main="residuals CD full model")
# hist(residuals(TL3), main="residuals CD full model")

#' ------------------------------------------------------------------------------------------------
#' Run sfa models
#' ------------------------------------------------------------------------------------------------

# Cobb Douglass
sfaCD1 <- sfa(log(yld) ~ log(N) + log(lab), data=db1)
sfaCD2 <- sfa(CD2_form, data=db1)

par(mfrow=c(1, 2))
hist(residuals(sfaCD1), main="residuals sfaCD core model")
hist(residuals(sfaCD2), main="residuals sfaCD full model")

# Translog
sfaTLenviron <- sfa(TL_form_environ, data=db1)
sfaTLGYGA <- sfa(TL_form_environ_GYGA, data=db1)

hist(residuals(sfaTLenviron), main="residuals sfaCD core model")
hist(residuals(sfaTLGYGA), main="residuals sfaCD full model")

# lrtest to determine what is a better fit


#' ------------------------------------------------------------------------------------------------
#' z variable analysis
#' ------------------------------------------------------------------------------------------------

sfaTL_z <- sfa(TL_form_fullz, data=db1)
hist(residuals(sfaTL_z), main="residuals sfaTL_z core model")
z_res <- summary(sfaTL_z)$mleParam
z_res <- z_res[grep("Z_", row.names(z_res)), ]
z_res <- as.data.frame(round(z_res, 3))


#' ------------------------------------------------------------------------------------------------
#' Endogeneity
#' ------------------------------------------------------------------------------------------------

db1$Nout <- ifelse(db1$noN == 1, 0, db1$logN)
db1$Nout <- ifelse(db1$Nout < 0, 0, db1$logN)

# estimate tobit model
# CHECK: Probably better to add relative price (WORKS!) (or prices for Pc and Pn separately) 
m <- tobit(Nout ~ dist_market + age + sex +
             ed_any + SOC + Pn +elevation +
             log(slope + 1) + rain_wq + phdum55_2_70 + phdum_gt70 + 
             extension, data=db1, left=0)
stargazer(m, type = "text")

# get residuals
rd <- resid(m, type = "deviance")

# have a look at fit of tobit model
plot(fitted(m), rd, main = "Fitted vs Residuals")
qqnorm(rd[!db1$Nout == 0]); abline(v=0, h=0); abline(v=2, h =2); abline(v=-2, h=-2)

# use the control function to test for endogeneity
# next step is get r into dataframe
indx <- m$na.action

db1 <- cbind(db1[-indx, ], rd)

# now run the translog again with only the environ vars
TL_form_environr <- paste(paste(deparse(formula(TL_form_environ), width.cutoff = 500),
                              "rd", sep=" + "), "-noN", sep=" ")

sfaTLenvironr <- sfa(TL_form_environr,
                     data = db1)

# run translog again with environ + GYGA variables
TL_form_environ_GYGAr <- paste(paste(deparse(formula(TL_form_environ_GYGA), width.cutoff = 500),
                                "rd", sep=" + "), "-noN", sep=" ")

sfaTLenvironGYGAr <- sfa(TL_form_environ_GYGAr,
                     data = db1)

# bootstrap SE - 500 runs ideally (takes 10 mins)
# not run except to get results for tables, otherwise
# takes too long!
# B <- 500
# coefM <- matrix(ncol=length(coef(sfa_TL_full_r)), nrow=B)
# system.time({
#   for(i in 1:B){
#     indx <- sample(1:nrow(db1.1), nrow(db1.1), replace=TRUE)
#     modl <- sfa(formula(sfa_TL_full_r), data=db1.1[indx, ])
#     coefM[i, ] <- coef(modl)
#   }
# })
# se <- apply(coefM, 2, sd)
# coef <- coef(sfa_TL_full_r)
# 
# N <- nrow(sfa_TL_full_r$dataTable)
# df <- N - length(coef)
# t <- coef/se
# 
# # conclusion from bootstrapping is that endogeneity
# # is present -> smaller N values
# p.value = round(2*pt(abs(t), df=df, lower=FALSE), 3)


#' ------------------------------------------------------------------------------------------------
#' table of results
#' ------------------------------------------------------------------------------------------------

# results from 4 models to compare
# 1. translog model without endogeneity or environmental variables

results <- full_join(sfa_table(sfaTLenviron, "Basic"),
                     sfa_table(sfaTLGYGA, "GYGA")) %>%
  full_join(., sfa_table(sfaTLenvironr, "Basic + r")) %>%
  full_join(., sfa_table(sfaTLenvironGYGAr, "GYGA + r"))
move <- which(results$parameter %in% c("gamma", "sigmaSq") )
results <- rbind(results[-move, ], results[move,])
results

#' ------------------------------------------------------------------------------------------------
#' z variable analysis - now with the residuals included
#' ------------------------------------------------------------------------------------------------




#' ------------------------------------------------------------------------------------------------
#' Find MPP and optimal Nitrogen. This cannot be done
#' analytically for a translog => numerical methods
#' ------------------------------------------------------------------------------------------------

# if there is no relative price we cannot
# work out the optimum
db1 <- db1[!is.na(db1$relprice),] # remove NA from relprice

#model that we want to find optimum for
TL_form_basic <- translog_form("log(yld)", translog_inputs)

# full specification including environmental inputs
TL_form_full <- paste(TL_form_basic,
                      paste(environmental, collapse=" + "),
                      paste(GYGA, collapse=" + "),
                      paste(dummies[!dummies %in% "noN"], collapse=" + "), "rd", sep=" + ") # leave noN variable out for now - uncleaer whether to include or not
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

par(mfrow=c(1, 1))

# nice example 
row <- db1.1[2, ]
N <- seq(1, 1000, by=0.1)
plot(N, Nopt_f(N, modl, row), type="l")
abline(h=0)

# not so nice example in which the function
# is never below zero -> no root!!
row <- db1.1[597, ]
N <- seq(-10, 1000, by=0.1)
plot(N, Nopt_f(N, modl, row), type="l", ylim=c(10, -100))
abline(h=0)

# have a look at fit to see if this can
# explain why sometimes there is no optimal solution
# clearly there are some large ish outliers
par(mfrow=c(1, 2))
plot(fitted(modl), residuals(modl))
plot(exp(modl$dataTable[, 3]), exp(fitted(modl)))
text(exp(modl$dataTable[, 3]), exp(fitted(modl)), 1:nrow(db1))

# function to iterate over. Finds optimal N 
# for each plot in the data. However, in order
# to use the uniroot function for this, we need
# to make sure that the lower and upper
# limit of the interval argument evaluate to
# values with opposite signs. Hence the while 
# and if statements in the following function

f <- function(i){
  print(i) # print iteration
  row <- db1.1[i, ] # select row
  
  # first find maximum (should be above zero in most cases)
  # this is used as the lower (positive) limit in the
  # uniroot interval
  res <- optimize(function(x) Nopt_f(x, modl, row), interval=c(0, 1000), maximum =TRUE)
  lower <- res$maximum
  if(res$objective < 0) return(NA)
  
  # choose a really large upper limit
  upper <- 20000
  ures <- uniroot(function(x) Nopt_f(x, modl, row), lower=lower, upper=upper)
  return(ures$root)
}

# iterate over each plot
db1.1$Npm <- sapply(1:nrow(db1.1), f)

# a few exploratory plots
par(mfrow=c(1, 1))
hist(db1.1$Npm, xlim=c(0, 700), breaks=500)
par(mfrow=c(1,2))
hist(db1.1$Npm[db1.1$noN==1], xlim=c(0, 500), breaks=500, main="no N")
hist(db1.1$Npm[db1.1$noN==0], xlim=c(0, 500), breaks=500, main="yes N")
mean(db1.1$Npm[db1.1$noN==1], na.rm=TRUE); mean(db1.1$Npm[db1.1$noN==0], na.rm=TRUE)

# CHECK 185 bad ones!
table(is.na(db1.1$Npm)) # 8 bad ones

db1.1$N <- Nobs # get the actual observed values back -> needed for MPP calculation
db1.1$Ndif <- db1.1$N - db1.1$Npm
hist(db1.1$Ndif, breaks=500, xlim=c(-700, 500))

# have a look at rows with no solution and
# see if they have anything in common
badrows <- db1.1[is.na(db1.1$Npm), ]

# fairly high relprice in badrows

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
  if(row$noN == 1){
    return(NA)
  } else {
  return(MPP_f(row))
  }
}

# get MPP values
db1.1$MPP <- sapply(1:nrow(db1.1), f2)

# histogram of MPP values
hist(db1.1$MPP, breaks=40)
mean(db1.1$MPP, na.rm=TRUE)

# plot against actual nitrogen use
# indicates that the higher users of
# nitrogen have lower MPP -> in accordance
# with theory
par(mfrow=c(1, 1))
with(db1.1, plot(N, MPP,  main="N vs. MPP"))


#' ------------------------------------------------------------------------------------------------
#' summarise at the zone level the MPP, optimal N
#' and actual N with number of observations
#' ------------------------------------------------------------------------------------------------

by_zone <- group_by(db1.1, ZONE) %>% summarise(
  n = n(),
  `N plots`  = n - sum(noN),
  N = round(mean(N[noN ==0], na.rm=TRUE),0),
  Nopt = round(mean(Npm, na.rm=TRUE), 0),
  MPP = round(mean(MPP, na.rm=TRUE), 0))


#' ------------------------------------------------------------------------------------------------
#' calculate the four types of yield
#' 1. Technical Efficiency Yield
#' 2. Economic Yield
#' 3. Feasible Yield
#' 4. Potential Yield
#' ------------------------------------------------------------------------------------------------

# frontier translog model
model <- modl

# 1. Technical efficiency yield is found using the
# output of the sfa model
# Observations where Npm cannot be calculated are removed
db3 <- db1.1 %>%
  rename(Y = yld) %>%
  mutate(
    Ycor = exp(as.numeric(fitted(model))+log(as.numeric(efficiencies(model)))), 
    err = Ycor-Y,
    TEY = exp(as.numeric(fitted(model))),
    TE = as.numeric(efficiencies(model)),
    resid = as.numeric(resid(model))
  ) %>%
  filter(!is.na(Npm))

table((db3$Y - db3$TEY) > 0) # number of plots where actual yield bigger than frontier yield

# 2. Economic yield is found by evaluating the frontier 
# function at the economically optimal nitrogen rate (Npm)

# introduce a cap on nitrogen
Npy <- 700 # check this

# Cap Npm
db3 <- mutate(db3, Npm = ifelse(Npm>Npy, Npy, Npm))

# evaluate frontier function at optimal level of
# N -> Npm, in order to get the economic yield
ysum <- sfaFormEval(modl)
ysum <- mgsub("log(N)", "log(Npm)", ysum)
db4 <- db3 %>% 
  mutate(EY = exp(with(., eval(parse(text=ysum)))))

# compare histograms of people who do
# and do not use nitrogen -> much the same
hist(db4$EY[db4$noN==1], col=rgb(1,0,0,0.5),
     main="EY of those using and not using N",
     xlab="N",
     breaks=40)
hist(db4$EY[db4$noN==0], col=rgb(0,0,1,0.5), add=T, breaks=40)


# 3. PFY: Feasible yield
# To improve this part, we could also argue that: (1) hybrid seeds are used, (2) pestices are used, (3) higher levels of capital and labour are used.
# We assume that all farmers use pesticides and increase assets and labour with 10%

# CHECK increase seed QUANTITY
# NEED TO SET THE IMPROVED SEED DUMMY TO 1 

# what should be multiplied by 1.1 here? loglab of log(lab*1.1)??
# Also note that N is evaluated at the cap of 120
library(qdap)
ysum <- sfaFormEval(modl)
ysum <- mgsub(dummies[-1],"1",ysum)
ysum <- mgsub(c("lab"), paste0(c("lab"), "*1.5"), ysum)
#ysum <- mgsub(c("log(lab)"), paste0("(", c("log(lab)"), "*1.5", ")"), ysum)

ysum <- mgsub("log(N)", "log(Npy)", ysum)

db5 <- db4 %>%
  mutate(PFY = exp(with(., eval(parse(text=ysum)))))

#check <- dplyr::select(db5, hhid, holder_id, parcel_id, field_id, surveyyear, ZONE, REGNAME, area, crop_count2, lat, lon, noN, yesN, loglab, lab, Npm, MPP, Y, PFY) 

# 4. PY: Potential yield
# Merge Yield potential with maize plot database
db6 <- db1 %>% dplyr::select(hhid, holder_id, parcel_id, field_id, PY = YW) %>% na.omit %>% 
  mutate(PY=PY*1000) %>% left_join(db5,.)


# A large number of plots have missing YW values because region is not covered by GYGA.
# We assume for the moment that country maximum water limited yield (Yw) is reasonable proxy for missing information.
# Might scale this down to see the effect.

GYGA_YW <- 18071.7857142857
db6 <- mutate(db6, PY = ifelse(is.na(PY), GYGA_YW, PY))

#####################################################
### Yield levels consistency check and correction ###
#####################################################

# Because of imputation of TY or measurement error, Yield (Y and Ycor),
# Technical efficiency yield (TEY), Economic yield (EY) and Feasible yield (UY) 
# can be higher than Potential yield (PYcor). We check for this.

X_Y_Ycor_check <- filter(db6, Ycor-Y<0)
X_PY_Y_check <- filter(db6, PY-Y<0)
X_PY_Y_cor_check <- filter(db6, PY-Ycor<0)
X_PY_TE_check <- filter(db6, PY-TEY<0)
X_PY_EY_check <- filter(db6, PY-EY<0)
X_PY_PFY_check <- filter(db6, PY-PFY<0)

# Compare different yield levels
# Picture shows that PFY is much to high for plots with the lowest PY.
# This is probably due to the uniform use of Npf of 120 N/ha.
# It would be better to have zone specific Npf values.
ggplot(data = db6, aes(y = PY, x = PFY)) +
  geom_point() +
  #geom_jitter(position=position_jitter(width=.1, height=0))+
  geom_abline(aes(Y = Ycor), slope=1, intercept=0) +
  coord_fixed() +
  scale_y_continuous(limits=c(0, 10000)) +
  scale_x_continuous(limits=c(0, 10000))

# Compare log(y) and resid
# Resid is defined as log(Y)-log(TEY) = error(v) - efficiency(u)
ggplot(data = db6, aes(y = log(Y), x = resid)) +
  geom_point() 

# Compare Sfa TA scores with mannually computed TEYG_s => identical as they should be
db6a <- mutate(db6, TEYG_s =Ycor/TEY)
ggplot(data = db6a, aes(y = TEYG_s, x = TE)) +
  geom_point()

#  We cap all values at PY because we consider this as an absolute potential and recalculate all gaps.
db7 <- mutate(db6, PFY = ifelse(PY-PFY<0, PY, PFY),
              EY = ifelse(PY-EY<0, PY, EY),
              TEY = ifelse(PY-TEY<0, PY, TEY),
              Ycor = ifelse(PY-Ycor<0, PY, Ycor),
              Y = ifelse(PY-Y<0, PY, Y))

#############################
### Yield gap calculation ###
#############################

# Calculate TYG using UY as reference
db8 <- db7 %>% 
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

# Consistency check of yield gaps.
# ERROR
X_ERROR_check <- filter(db8, ERROR_l<0) # Half of observation has a negative error which is what would be expected
mean(db8$ERROR_l)
mean(db8$ERROR_s)

# TEYG
X_TEYG_check <- filter(db8, TEYG_l<0) # Should be zero
mean(db8$TEYG_s)

# EYG
# A number of plots will have to decrease N use Npm < N. In several cases also plots that do no use N
# will have lower Y when they start using N. This is because there yield can be located above the frontier (based on fertilizer users) because of the positive effect of noN.
# If we believe that these plots are structurally different and do not use fertilizer because of better soils, they will in fact use too much N and have to decrease.
X_EYG_check <- filter(db8, EYG_l<0)        
mean(db8$EYG_s)

# EUYG
# A number of plots have negative EUYG_l because Npm is larger than Npy, the nitrogen that is required to achieve Potential yield (Yw).
# We have corrected this so check should be 0.
# CHECK: we stil find negative values. It is possible that because of the interaction with labour >N results in <y? TO BE CHECKED
X_EUYG_check <- filter(db8, EUYG_l<0)        
mean(db8$EUYG_s)
check <- dplyr::select(X_EUYG_check, hhid, holder_id, parcel_id, field_id, surveyyear, ZONE, REGNAME, area, crop_count2, lat, lon, noN, yesN, loglab, lab, Npm, MPP, Y, 
                       PFY, EY, EUYG_l) 

# TYG
X_TYG_check <- filter(db8, TYG_l<0)        
mean(db8$TYG_s)

#YG
X_YG_check <- filter(db8, YG_l<0)        
YG_check2 <- filter(db8, YG_l_Ycor<0)        

# Check if separate yield gaps add up to total yield gap
Overall_check <- db8 %>%
  mutate(check_l = YG_l/(ERROR_l + TEYG_l + EYG_l + EUYG_l + TYG_l), # Note that for a small number of observatios YG_l=0 resulting in 0/0 which is NaN
         check_s = YG_s/(ERROR_s * TEYG_s * EYG_s * EUYG_s * TYG_s),
         check_l2 = YG_l_Ycor/(TEYG_l + EYG_l + EUYG_l + TYG_l),
         check_s2 = YG_s_Ycor/(TEYG_s * EYG_s * EUYG_s * TYG_s))
summary(Overall_check)


# Create database with relevant variables for further analysis
db9 <- dplyr::select(db8, hhid, holder_id, parcel_id, field_id, ZONE, REGNAME, surveyyear, lat, lon, crop_count2, area, Npm, yesN, Y, N, Ycor, TEY, EY, PFY, PY, ERROR_l, ERROR_s, TEYG_l, TEYG_s, EYG_l, EYG_s, 
                     EUYG_l, EUYG_s, TYG_l, TYG_s, YG_l, YG_s, YG_l_Ycor, YG_s_Ycor)

saveRDS(db9, "Cache/db9.rds")
saveRDS(results, "Cache/results_table.rds")
saveRDS(dbsum, "Cache/dbsum.rds")
saveRDS(by_zone, "Cache/by_zone.rds")
saveRDS(z_res, "Cache/z_res.rds")
