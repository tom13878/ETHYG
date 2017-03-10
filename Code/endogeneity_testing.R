#'==============================================================================
#' Project:  IMAGINE ETH
#' Subject:  Endogeneity tests for CD and TL functions
#' Author:   Michiel van Dijk & Tom Morley
#' Contact:  michiel.vandijk@wur.nl, Tomas.morley@wur.nl
#' Output:   Results table of endogeneity tests
#'==============================================================================

# get packages
library(pacman)
p_load(char=c("rprojroot", "boot", "frontier", "dplyr", "AER"), install=TRUE)

# get root path
root <- find_root(is_rstudio_project)

# read in data to run models on
db1 <- readRDS(file.path(root, "Cache/db1.rds"))
db1 <- unique(db1)
db1 <- filter(db1, lab < 2000, logN >= 0)

# make some new variables for convenience.
db1$logNsq <- (db1$logN)^2
db1$loglabsq <- db1$loglab^2
db1$logNloglab <- db1$logN * db1$loglab
db1$logslope <- log(db1$slope + 1)
db1$agesq <- db1$age^2

# -------------------------------------
# translog production function with 
# endogeneity

# there are some missing values so we exclude these from the data now to avoid
# problems later
db1 <- as.data.frame(model.matrix(~ -1 + logyld + yesN + logN + loglab + loglabsq +
                      logNloglab + Pn + dist_market + sex + age + agesq +
                      ed_any + logslope + elevation +
                      logarea + crop_count2, data=db1))

# for comparison we the following results we begin
# with a simple translog function as a benchmark
TL_simple <- lm(logyld ~ logN + loglab + loglabsq +
                     logNloglab + yesN + logslope + elevation +
                     logarea + crop_count2, data=db1)


# first stage linear model
stage1_lm <- lm(logN ~ Pn + dist_market + sex + age + agesq + ed_any +
                       loglab + loglabsq + logslope + elevation +
                       logarea + crop_count2, data=db1)

# get residuals
db1$v <- rstandard(stage1_lm) 

# second stage TL model with residuals
stage2_TL_lm <- lm(logyld ~ logN + loglab + loglabsq +
                        logNloglab + yesN + logslope + elevation +
                        logarea + crop_count2 + v, data=db1)

# bootstrap results - need to include ALL stages in bootstrap
refit_TL_lm <- function(data, indx){
  dat <- data[indx,]
  fs <- lm(formula(stage1_lm), data=dat)
  dat$v <- rstandard(fs) 
  coef(lm(formula(stage2_TL_lm), data=dat))
}

TL_lm_boot <- boot(db1, refit_TL_lm, R = 500)

# first stage tobit model
stage1_tob <- tobit(logN ~ Pn + dist_market + sex + age + agesq + ed_any +
                       loglab + loglabsq + logslope + elevation +
                       logarea + crop_count2, data=db1)

# calculate the generalized residual (Greene)
d1 <- 1 - db1$yesN
d2 <- db1$yesN
sigma <- stage1_tob$scale
theta <- 1/sigma
mills <- -dnorm(-fitted(stage1_tob))/pnorm(-fitted(stage1_tob))
db1$v <- d1 * mills + d2 * (theta * db1$logN - fitted(stage1_tob))

# put v into the second stage TL model
stage2_TL_tob <- lm(logyld ~ logN + loglab + loglabsq +
                     logNloglab + logslope + elevation +
                     logarea + crop_count2 + v, data=db1)

# bootstrap results: note that for some reason we cannot use the
# formula(modl) trick with the tobit function -> have to write out everything
# whole
refit_TL_tob <- function(data, indx){
  dat <- data[indx, ]
  fs <- tobit(logN ~ Pn + dist_market + sex + age + agesq + ed_any +
                loglab + loglabsq + logslope + elevation +
                logarea + crop_count2, data=dat)
  d1 <- 1 - dat$yesN
  d2 <- dat$yesN
  sigma <- fs$scale
  theta <- 1/sigma
  mills <- -dnorm(-fitted(fs))/pnorm(-fitted(fs))
  dat$v <- d1 * mills + d2 * (theta * dat$logN - fitted(fs))
  coef(lm(formula(stage2_TL_tob), data=dat))
}

TL_tob_boot <- boot(db1, refit_TL_tob, R = 500)

# tobit results are not comparable to those from
# an OLS model unless we first calcualte the APE
# as follows. Wooldridge (2010) pp. 674

sigma <- stage1_tob$scale
X <- model.matrix(stage1_tob)
n <- nrow(X)
scale_factor <- (1/n) * sum(pnorm(((X %*% coef(stage1_tob)))/sigma))
APE <- coef(stage1_tob) * scale_factor

# however, to get SEs for the APE we also need to
# bootstrap the making of the APE!!!
refit_APE <- function(data, indx){
  dat <- data[indx, ]
  fs <- tobit(logN ~ Pn + dist_market + sex + age + agesq + ed_any +
                loglab + loglabsq + logslope + elevation +
                logarea + crop_count2, data=dat)
  sigma <- fs$scale
  X <- model.matrix(fs)
  n <- nrow(X)
  scale_factor <- (1/n) * sum(pnorm(((X %*% coef(fs)))/sigma))
  APE <- coef(fs) * scale_factor
  APE
}

# results from bootstrap
APE_boot <- boot(db1, refit_APE, R = 500)

# -------------------------------------
# construct a table with the results

# first stage table
OLS <- summary(stage1_lm)$coef
Tobit <- summary(stage1_tob)$coef
APE <- summary(APE_boot)
FS_tab <- cbind(OLS[, 1:2], Tobit[-nrow(Tobit), 1:2], APE[, c(2, 4)])
FS_tab <- round(FS_tab, 3)
names(FS_tab) <- c("OLS", "SE", "Tobit", "SE", "APE", "BSE")

# second stage table
TL_simple <- rbind(summary(TL_simple)$coef[-7, ], v=NA)
TL_lm <- summary(stage2_TL_lm)$coef[-7, ] # remove yesN because it is not in both models
TL_tob <- summary(stage2_TL_tob)$coef

# bootstrapped SEs
TL_lm_SE <- summary(TL_lm_boot)[-7, ]
TL_tob_SE <- summary(TL_tob_boot)
SS_tab <- as.data.frame(cbind(TL_simple[, 1:2], TL_lm[, 1:2], TL_lm_SE[, 4], TL_tob[, 1:2], TL_tob_SE[, 4]))
SS_tab <- round(SS_tab, 3)
names(SS_tab) <- c("TL", "SE", "TL-LM", "SE", "BSE", "TL-Tob", "SE", "BSE")

# save tables for later use
saveRDS(FS_tab, file.path(root, "report/tables/FS_endog_tab.rds"))
saveRDS(SS_tab, file.path(root, "report/tables/SS_endog_tab.rds"))

# finally, consider a heckit sample selection model
# with endogeneity. remember we are asking what is the
# effect of nitrogen per hectare on maize yield.

# step 1 is to do a probit model and get mills ratio
probit <- glm(yesN ~ Pn + dist_market + sex + age + agesq + ed_any +
      loglab + loglabsq + logslope + elevation +
      logarea + crop_count2, data=db1, family=binomial("probit"))
summary(probit)
db1$mills <- dnorm(fitted(probit))/pnorm(fitted(probit)) 

# step 2 make a sub sample selection based on seeing N or not
sub_db1 <- db1[db1$yesN == 1, ]

# step 3 estimate the new equation by 2sls or control function
stage1_lm <- lm(logN ~ Pn + dist_market + sex + age + agesq + ed_any +
                  loglab + loglabsq + logslope + elevation +
                  logarea + crop_count2 + mills, data=sub_db1)

sub_db1$v <- rstandard(stage1_lm)

stage2_TL_lm <- lm(logyld ~ logN + loglab + loglabsq +
                     logNloglab  + logslope + elevation +
                     logarea + crop_count2 + v + mills, data=sub_db1)

# and bootstrap SE
refit_TL_lm <- function(data, indx){
  dat <- data[indx,]
  fs <- lm(formula(stage1_lm), data=dat)
  dat$v <- rstandard(fs) 
  coef(lm(formula(stage2_TL_lm), data=dat))
}

TL_lm_boot <- boot(sub_db1, refit_TL_lm, R = 500)