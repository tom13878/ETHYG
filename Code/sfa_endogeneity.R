#'==============================================================================
#' Project:  IMAGINE ETH
#' Subject:  sfa with endogeneity
#' Author:   Michiel van Dijk & Tom Morley
#' Contact:  michiel.vandijk@wur.nl, Tomas.morley@wur.nl
#' Output:   Results table of sfa endogeneity 
#'==============================================================================

# get packages
library(pacman)
p_load(char=c("rprojroot", "ucminf", "optimx", "dplyr", "frontier", "boot"), install=TRUE)

# get root path
root <- find_root(is_rstudio_project)

# source in sfa functions
source(file.path(root, "Code/sfa_functions.R"))

# read in data to run models on
db1 <- readRDS(file.path(root, "Cache/db1.rds"))

# -------------------------------------
# translog production function with 
# endogeneity

# there are some missing values so we exclude these from the data now to avoid
# problems later
db1 <- as.data.frame(model.matrix(~ -1 + logyld +  logN + loglab  + logseed +
                                    logNsq + loglabsq + logseedsq +
                                    logNlab + logNseed + loglabseed +
                                    Pn + dist_market + cost2large_town,
                                  data=db1))


# step 1: test if endogeneity exists
stage1 <- lm(logN ~ Pn + cost2large_town + dist_market, data=db1)
v <- rstandard(stage1)
stage2 <- lm(logyld ~ logN + loglab + logseed +
             logNsq + loglabsq + logseedsq +
             logNlab + logNseed + loglabseed + v, data=db1)

# bootstrap SE
refit <- function(data, indx){
  dat <- data[indx,]
  fs <- lm(formula(stage1), data=dat)
  v <- residuals(fs) 
  coef(lm(formula(stage2), data=dat))
}

# bootstrap results
boot_out <- summary(boot(db1, refit, R = 500))

# p.values
N <- nrow(db1)
df <- N - length(stage2$coef)
t <- stage2$coef/boot_out[, 4]
pval <- round(2*pt(abs(t), df=df, lower=FALSE), 3)

# and for comparison just an ordinary ols translog
olstl <- lm(logyld ~ logN + loglab + logseed +
                        logNsq + loglabsq + logseedsq +
                        logNlab + logNseed + loglabseed, data=db1)


# Step 2: liml estimation using likelihood from APS 2016

# basic frontier model for comparison
sf <- summary(sfa(logyld ~ logN + loglab + logseed +
                    logNsq + loglabsq + logseedsq +
                    logNlab + logNseed + loglabseed, data=db1))

# data
# outcome of structural equation
Y <- db1$logyld

# endogenous variable(s)
X2 <- model.matrix( ~ -1 + logN, data=db1)

# variables we want in the structural equation but not 
# the reduced equation
X1 <- model.matrix(~ -1 + loglab + logseed +
                     logNsq + loglabsq + logseedsq +
                     logNlab + logNseed + loglabseed,
                   data = db1)

# full X ( all structural variables )
X <- cbind(X2, X1)

# instruments/variables only in reduced equation
W <- model.matrix(~ -1 + Pn + cost2large_town + dist_market,
                  data = db1)

# complete set of regressors for reduced equation
Z <- cbind(1, W)
N <- length(Y)

# initial parameters for optimisation - note
# that parameter sigma2v must be greater than
# Sigmavn * (1/Sigmann) * Sigmavn, otherwise 
# sigma2c will be negative and this is not
# allowed.
n <- (ncol(X) + ncol(Z) + 5)
pars <- rep(0.2, n)
pars[ncol(X) + 2] <- 2

# test out the function based on
# starting parameters - starting parameters
# must give a non NAN answer
liml2(pars, X, X2, Y, Z, N)

# optimx function
stats <- optimx(pars, liml2, NULL,
                hess=NULL,
                lower=-Inf,
                upper=Inf,
                method="BFGS",
                itnmax=NULL,
                hessian=FALSE,
                control = list(trace=2, maxit=200),
                X, X2, Y, Z, N)

# get parameters (same order as they appear in the liml1 function)

# output from optimx
alpha <- stats[1]
beta <- stats[2:(1 + ncol(X))]
sigma2v <- stats[(2 + ncol(X))]
sigma2u <- stats[(3 + ncol(X))]
pos <- (4+ncol(X))
Sigmavn <- stats[pos]
Sigmann <- stats[pos + 1]
Pi <- stats[(pos + 2): (pos + 1 + ncol(Z))]

# match up beta parameters with the names from
# the X matrix (second stage)
unlist(c(alpha, beta)) # parameters from liml second stage
dimnames(X)[[2]] # order in which parameters appear from liml output
sf$mleParam[, 1] # normal (no endogeneity) SF ML estimates
sf$olsParam[, 1] # ols estimates
c(sigma2v, sigma2u) # variance and relative variance estimates

# match up Pi parameters with the names from the Z
# matrix (first stage)
Pi # parameters from liml first stage
dimnames(Z)[[2]] # order in which parameters appear from liml output
stage1 # ols normal reduced model

# ==============================================================================
# make tables of results comparing estimations.

# testing endogeneity tables
# stage 1 - ols vs liml ML
stage1_tab <- round(as.data.frame(summary(stage1)$coef[, -4]), 3)
liml1_tab <- as.numeric(round(Pi, 3))
stage1_tab <- data.frame(stage1_tab, liml=liml1_tab)
saveRDS(stage1_tab, file.path(root, "Cache/endog1.rds"))

# stage 2 - control function vs std ols and test for endog
stage2_tab <- round(as.data.frame(summary(stage2)$coef[, -4]), 3)
stage2_tab$`Std. Error` <- boot_out$bootSE
stage2_tab$`t value` <- t
olstl_tab <- round(as.data.frame(summary(olstl)$coef[, -4]), 3)
olstl_tab <- rbind(olstl_tab, v=c(NA, NA, NA))
stage2_tab <- cbind(olstl_tab, stage2_tab)
names(stage2_tab)[4:6] <- c("CF Est", "CF SE", "CF p.val")
saveRDS(stage2_tab, file.path(root, "Cache/endog2.rds"))

# accounting for endogeneity tables
sf_tab <- round(as.data.frame(sf$mleParam), 3)

# compare normal sf with the liml estimates
sigma2 <- sigma2u + sigma2v
gamma <- sigma2u/sigma2
liml_tab <- as.numeric(round(unlist(c(alpha, beta, sigma2, gamma)),3))
sf_tab <- data.frame(sf_tab, LIML=liml_tab)
saveRDS(sf_tab, file.path(root, "Cache/endog3.rds"))
