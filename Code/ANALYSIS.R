#'==============================================================================
#' Project:  IMAGINE ETH
#' Subject:  Analysis file
#' Author:   Michiel van Dijk & Tom Morley
#' Contact:  michiel.vandijk@wur.nl, Tomas.morley@wur.nl
#' Output:   yield gap breakdown 
#'==============================================================================

# get packages
library(pacman)
p_load(char=c("dplyr", "rprojroot", "frontier"), install=TRUE)
root <- find_root(is_rstudio_project)

# get data
db1 <- readRDS(file.path(root, "Cache/db1.rds"))

# there is no predict function in the sfa package
# but we would like one - make our own and source
# it in
source(file.path(root, "Code/predict.sfa.R"))

# run the model
sf11x9 <- sfa(logyld ~ logN + loglab + logseed +
                logNsq + loglabsq + logseedsq +
                logN:loglab + logN:logseed +
                loglab:logseed + logarea + phdum55_2_70 +
                crop_count2 + dumoxen + SOC2 + logslope +
                elevation + GGD + AI + TS|
                -1 + age + sex + ed_any + title +
                extension + credit + dist_market +
                popEA + logarea_tot, data=db1)

# we want to evaluate the frontier function
# and avoid the Z variables
zidx <- grep("Z_", names(coef(sf11x9)))
xvars <- names(coef(sf11x9))[1:(min(zidx)-1)]
X <- sf11x9$dataTable[, xvars]
X <- as.data.frame(X)
xcoef <- coef(sf11x9)[1:(min(zidx)-1)]
relprices <- db1$relprice[sf11x9$validObs]

# function to calculate the MPP
calc_mpp <- function(N, row){
  row$logN <- log(N) 
  logY <- as.matrix(row) %*% xcoef
  Y <- exp(logY)
  MPP <- with(row, ((xcoef["logN"] + 
                       2*xcoef["logNsq"]*logN +
                       xcoef["logN:loglab"]*loglab +
                       xcoef["logN:logseed"]*logseed )*
                      (Y/N)))
  MPP
}

# function to calculate the mpp less
# the relative price 
mpp_price <- function(N, row, relprice){
  calc_mpp(N, row) - relprice
}

# function to calculate the economically
# optimal level of nitrogen
econ_opt <- function(i){
  
  # get a plot observation and relative
  # price
  row <- X[i, ]
  relprice <- relprices[i]
  
  # first we need a point above zero
  # but before the root to act as the
  # lower limit of the root finding 
  # function - We cap at 700 because
  # even on the most productive farms
  # 700 kg/ha of nitrogen is excessive
  lower <- tryCatch(optimize(function(x) mpp_price(x, row, relprice),
                       interval=c(1, 700),
                       maximum =TRUE)$maximum, error=function(err) NA)
  
  # then we need to find the root
  # i.e. the point at which the
  # nitrogen level crosses zero. 
  # this is the economically
  # constrained level of nitrogen
  tryCatch(uniroot(function(x) {mpp_price(x, row, relprice)},
                   interval=c(lower, 700))$root,
           error=function(err) NA)
}

# calculate the economically
# optimal level of nitrogen per plot
db2 <- db1[sf11x9$validObs, ]
db2$Npm <- sapply(1:nrow(X), econ_opt)

# we can calculate a mpp wrt nitrogen
# at the observed nitrogen values.
# For values of N that are less than 1
# the MPP will be negative and large
# because the formula for the mpp of a
# translog involves the term 1/N
# for a value of N = 0.01, for example,
# this will result in an multiplicative
# of 100.
# of nitrogen that are 0, this will be 
# negative. DOn;t forget the -1 we added
# to N
db2$mpp <- calc_mpp((exp(X$logN) - 1) , X)

model <- sf11x9

# 1. Technical efficiency yield is found using the
# output of the sfa model
# Observations where Npm cannot be calculated are removed
db2 <- db2 %>%
  rename(Y = yld) %>%
  mutate(
    Ycor = exp(as.numeric(fitted(model))+log(as.numeric(efficiencies(model)))), 
    err = Ycor-Y,
    TEY = exp(as.numeric(fitted(model))),
    TE = as.numeric(efficiencies(model)),
    resid = as.numeric(resid(model))
)

# 2. Economic yield is found by evaluating the frontier 
# function at the economically optimal nitrogen rate (Npm)
# This means we need to swap out the N (logN) variable
# for the Npm variable, BUT we also need a way of doing
# this for the interaction terms that also involve N
model_vars <- names(sf11x9$olsParam)[-c(1, length(names(sf11x9$olsParam)))]
model_vars <- model_vars[-grep(":", model_vars)] # get rid of interaction terms
predict_dat <- db2[, model_vars]
predict_dat <- mutate(predict_dat,
                      logN = log(db2$Npm),
                      logNsq = logN^2)

# now make the prediction predict.sfa
# function is not made to handle NA values
# so to keep order and compare with other yield
# measures we probably want to set NA values to
# zero temporarily 
predict_dat$logN[is.na(predict_dat$logN)] <- 0
predict_dat$logNsq[is.na(predict_dat$logNsq)] <- 0
db2$EY <- exp(predict.sfa(sf11x9, predict_dat))
db2$EY[predict_dat$logN == 0] <- NA

# 3. PFY: Feasible yield
# evaluate frontier function at N = 150
# Increase labour and seed rate by 10%
# turn on all dummies.
predict_dat2 <- mutate(predict_dat,
                       logN = log(400),
                       logNsq = logN^2,
                       loglab = loglab + log(1.5),
                       loglabsq = loglab^2,
                       logseed = logseed + log(1.5),
                       logseedsq = logseed^2,
                       dumoxen = 1)

# make prediction
db2$PFY <- exp(predict.sfa(sf11x9, predict_dat2))

# 4. Potential yield
db2$PY <- db2$YW * 1000

# A large number of plots have missing YW values because region
# is not covered by GYGA. We assume for the moment that country
# maximum water limited yield (Yw) is a reasonable proxy for
# missing information. Might scale this down to see the effect.

GYGA_YW <- max(db2$PY, na.rm=TRUE)
db2 <- mutate(db2, PY = ifelse(is.na(PY), GYGA_YW, PY))

# join all the yield measures into a single
# dataframe
#  We cap all values at PY because we consider this as an absolute potential and recalculate all gaps.
db2 <- mutate(db2, PFY = ifelse(PY-PFY<0, PY, PFY),
              EY = ifelse(PY-EY<0, PY, EY),
              TEY = ifelse(PY-TEY<0, PY, TEY),
              Ycor = ifelse(PY-Ycor<0, PY, Ycor),
              Y = ifelse(PY-Y<0, PY, Y))

# Calculate TYG using UY as reference
db2 <- db2 %>% 
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
X_ERROR_check <- filter(db2, ERROR_l<0) # Half of observation has a negative error which is what would be expected
mean(db2$ERROR_l, na.rm=TRUE)
mean(db2$ERROR_s, na.rm=TRUE)

# TEYG
X_TEYG_check <- filter(db2, TEYG_l<0) # Should be zero
mean(db2$TEYG_s, na.rm=TRUE)

# EYG
# A number of plots will have to decrease N use Npm < N. In several cases also plots that do no use N
# will have lower Y when they start using N. This is because there yield can be located above the frontier (based on fertilizer users) because of the positive effect of noN.
# If we believe that these plots are structurally different and do not use fertilizer because of better soils, they will in fact use too much N and have to decrease.
X_EYG_check <- filter(db2, EYG_l<0)        
mean(db2$EYG_s)

# EUYG
# A number of plots have negative EUYG_l because Npm is
# larger than Npy, the nitrogen that is required to
# achieve Potential yield (Yw). We have corrected this
# so check should be 0.
# CHECK: we stil find negative values. It is possible
# that because of the interaction with labour >N results
# in <y? TO BE CHECKED
X_EUYG_check <- filter(db2, EUYG_l<0)        
mean(db2$EUYG_s)
check <- select(X_EUYG_check, hhid, holder_id, parcel_id, field_id,
                surveyyear, ZONE, REGNAME, area, crop_count2,
                lat, lon, noN, yesN, loglab, lab, Npm, mpp, Y, 
                       PFY, EY, EUYG_l) 

# TYG
X_TYG_check <- filter(db2, TYG_l<0)        
mean(db2$TYG_s, na.rm=TRUE)

#YG
X_YG_check <- filter(db2, YG_l<0)        
YG_check2 <- filter(db2, YG_l_Ycor<0)        

# Check if separate yield gaps add up to total yield gap
Overall_check <- db2 %>%
  mutate(check_l = YG_l/(ERROR_l + TEYG_l + EYG_l + EUYG_l + TYG_l), # Note that for a small number of observatios YG_l=0 resulting in 0/0 which is NaN
         check_s = YG_s/(ERROR_s * TEYG_s * EYG_s * EUYG_s * TYG_s),
         check_l2 = YG_l_Ycor/(TEYG_l + EYG_l + EUYG_l + TYG_l),
         check_s2 = YG_s_Ycor/(TEYG_s * EYG_s * EUYG_s * TYG_s))
summary(Overall_check)


# Create database with relevant variables for further analysis
db3 <- select(db2, hhid, holder_id, parcel_id, field_id, ZONE,
              REGNAME, surveyyear, lat, lon, crop_count2, area,
              Npm, yesN, Y, N, Ycor, TEY, EY, PFY, PY, ERROR_l,
              ERROR_s, TEYG_l, TEYG_s, EYG_l, EYG_s, EUYG_l,
              EUYG_s, TYG_l, TYG_s, YG_l, YG_s, YG_l_Ycor, YG_s_Ycor)

# save db3 for further analysis
saveRDS(db3, "Cache/db3.rds")
