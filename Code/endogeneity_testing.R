#'==============================================================================
#' Project:  IMAGINE ETH
#' Subject:  Endogeneity tests for TL functions
#' Author:   Michiel van Dijk & Tom Morley
#' Contact:  michiel.vandijk@wur.nl, Tomas.morley@wur.nl
#' Output:   Results table of endogeneity tests
#'==============================================================================

# set project root
library(rprojroot)
root <- find_root(is_rstudio_project)

# load packages
library(pacman)
p_load(frontier)

# get data
db1 <- readRDS(file.path(root, "Cache/db1.rds"))
db1 <- unique(db1)

# get square and interaction terms for translog
db1$logNsq <- db1$logN^2
db1$loglabsq <- db1$loglab^2
db1$logseedsq <- db1$logseed^2
db1$logarea_tot <- log(db1$area_tot)
db1$agesq <- db1$age^2

# Also divide GGD and AI by 1000
# to get a more interpretable coef
db1$GGD <- db1$GGD/1000
db1$AI <- db1$AI/1000

# -------------------------------------
# translog production function with 
# endogeneity

# first stage linear model
stage1_lm <- lm(logN ~ relprice + sex + age + ed_any +
                  logarea_tot + credit + extension +
                  title + elevation + logslope +
                  logarea + SOC2, data=db1)

# make a new database which has only those
# terms that we have complete cases for
db2 <- db1[-stage1_lm$na.action, ]

# get residuals
db2$v <- residuals(stage1_lm) 

# try doing stepwise selection
stage2_TL_lm <- lm(logyld ~ logN + loglab + logseed +
                logNsq + loglabsq + logseedsq +
                logN:loglab + logN:logseed +
                loglab:logseed + logarea + phdum55_2_70 +
                crop_count2 + dumoxen + SOC2 + logslope +
                elevation + GGD + AI + TS + v, data=db2)

summary(stage2_TL_lm)

# bootstrap results - need to include ALL stages in bootstrap
refit_TL_lm <- function(data, indx){
  dat <- data[indx,]
  fs <- lm(formula(stage1_lm), data=dat)
  dat$v <- residuals(fs) 
  coef(lm(formula(stage2_TL_lm), data=dat))
}

TL_lm_boot <- boot(db2, refit_TL_lm, R = 500)
TL_lm_SE <- summary(TL_lm_boot)

# we need to calcualte p.values using the BSEs
# lm first stage
N <- nrow(db2)
df <- N - length(stage2_TL_lm$coef)
t <- stage2_TL_lm$coef/TL_lm_SE[, 4]
p.value_lm <- round(2*pt(abs(t), df=df, lower=FALSE), 3)


# first stage tobit model
stage1_tob <- tobit(logN ~ relprice + sex + age + ed_any +
                      logarea_tot + credit + extension +
                      title + elevation + logslope +
                      logarea + SOC2, data=db2)

# calculate the generalized residual (Greene)
d1 <- 1 - db2$yesN
d2 <- db2$yesN
sigma <- stage1_tob$scale
theta <- 1/sigma
mills <- -dnorm(-fitted(stage1_tob))/pnorm(-fitted(stage1_tob))
db2$v <- d1 * mills + d2 * (theta * db2$logN - fitted(stage1_tob))

# put v into the second stage TL model
stage2_TL_tob <- lm(formula(stage2_TL_lm), data=db2)

# bootstrap results: note that for some reason we cannot use the
# formula(modl) trick with the tobit function -> have to write out everything
# whole
refit_TL_tob <- function(data, indx){
  dat <- data[indx, ]
  fs <- tobit(logN ~ relprice + sex + age + ed_any +
                logarea_tot + credit + extension +
                title + elevation + logslope +
                logarea + SOC2, data=dat)
  d1 <- 1 - dat$yesN
  d2 <- dat$yesN
  sigma <- fs$scale
  theta <- 1/sigma
  mills <- -dnorm(-fitted(fs))/pnorm(-fitted(fs))
  dat$v <- d1 * mills + d2 * (theta * dat$logN - fitted(fs))
  coef(lm(formula(stage2_TL_tob), data=dat))
}

TL_tob_boot <- boot(db2, refit_TL_tob, R = 500)
TL_tob_SE <- summary(TL_tob_boot)

# we need to calcualte p.values using the BSEs
# lm first stage
N <- nrow(db2)
df <- N - length(stage2_TL_tob$coef)
t <- stage2_TL_tob$coef/TL_tob_SE[, 4]
p.value_tob <- round(2*pt(abs(t), df=df, lower=FALSE), 3)

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
  fs <- tobit(logN ~ relprice + sex + age + ed_any +
                logarea_tot + credit + extension +
                title + elevation + logslope +
                logarea + SOC2, data=dat)
  sigma <- fs$scale
  X <- model.matrix(fs)
  n <- nrow(X)
  scale_factor <- (1/n) * sum(pnorm(((X %*% coef(fs)))/sigma))
  APE <- coef(fs) * scale_factor
  APE
}

# results from bootstrap
APE_boot <- boot(db2, refit_APE, R = 500)
APE_SE <- summary(APE_boot)

# work out p-values
N <- nrow(db2)
df <- N - length(APE)
t <- APE/APE_SE[, 4]
p.value_APE <- round(2*pt(abs(t), df=df, lower=FALSE), 3)

# -------------------------------------
# construct a table with the results

# first stage table
OLS <- summary(stage1_lm)$coef
Tobit <- summary(stage1_tob)$coef
APE <- summary(APE_boot)
FS_tab <- cbind(OLS[, 1:2], Tobit[-nrow(Tobit), 1:2], APE[, c(2, 4)])
FS_tab <- round(FS_tab, 3)
names(FS_tab) <- c("OLS", "SE", "Tobit", "SE", "APE", "BSE")

# add stars to tables
# ols first stage
stars <- cut(OLS[, 4], breaks=c(0, 0.01, 0.05, 1),
             labels=c("**", "*", ""), include.lowest = TRUE)
FS_tab$OLS <- paste0(FS_tab$OLS, stars)

# tobit first stage
stars <- cut(Tobit[-nrow(Tobit), 4], breaks=c(0, 0.01, 0.05, 1),
             labels=c("**", "*", ""), include.lowest = TRUE)
FS_tab$Tobit <- paste0(FS_tab$Tobit, stars)


# we can compare the lm and tobit first stage models to see which
# fits better in terms of the R squared. However, wooldridge (2010)
# pp. 680 points out that the OLS estimates are maximising the 
# R2 whereas the tobit estimates are maximising a loglikelihood
# however this still gives us an idea of which one is preferred.
# tobit first stage R-squared equivalent
yi <- pnorm(fitted(stage1_tob)/stage1_tob$scale)*fitted(stage1_tob) +
  stage1_tob$scale * dnorm(fitted(stage1_tob)/stage1_tob$scale)
Rsq_tob <- cor(yi, db2$logN)^2

Rsq <- round(c(summary(stage1_lm)$adj.r.squared, NA, Rsq_tob, NA, NA, NA), 3)
FS_tab <- rbind(FS_tab, Rsq)
row.names(FS_tab)[nrow(FS_tab)] <- "R-squared"

# second stage table
TL_lm <- summary(stage2_TL_lm)$coef
TL_tob <- summary(stage2_TL_tob)$coef

# bootstrapped SEs
TL_lm_SE <- summary(TL_lm_boot)
TL_tob_SE <- summary(TL_tob_boot)
SS_tab <- as.data.frame(cbind(TL_lm[, 1:2], TL_lm_SE[, 4],
                              TL_tob[, 1:2], TL_tob_SE[, 4]))
SS_tab <- round(SS_tab, 3)
names(SS_tab) <- c("TL-LM", "SE", "BSE", "TL-Tob", "SE", "BSE")

# add stars to tables
# lm -> translog
stars <- cut(p.value_lm, breaks=c(0, 0.01, 0.05, 1),
             labels=c("**", "*", ""), include.lowest = TRUE)
SS_tab$`TL-LM` <- paste0(SS_tab$`TL-LM`, stars)

# tobit -> translog
stars <- cut(p.value_tob, breaks=c(0, 0.01, 0.05, 1),
             labels=c("**", "*", ""), include.lowest = TRUE)
SS_tab$`TL-Tob` <- paste0(SS_tab$`TL-Tob`, stars)

# save tables for later use
saveRDS(FS_tab, file.path(root, "report/tables/FS_endog_tab.rds"))
saveRDS(SS_tab, file.path(root, "report/tables/SS_endog_tab.rds"))
