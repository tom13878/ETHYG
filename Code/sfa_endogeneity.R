#'==============================================================================
#' Project:  IMAGINE ETH
#' Subject:  sfa with endogeneity
#' Author:   Michiel van Dijk & Tom Morley
#' Contact:  michiel.vandijk@wur.nl, Tomas.morley@wur.nl
#' Output:   Results table of sfa endogeneity 
#'==============================================================================

# get packages
library(pacman)
p_load(char=c("rprojroot", "optimx", "dplyr", "frontier"), install=TRUE)

# get root path
root <- find_root(is_rstudio_project)

# source in sfa functions
source(file.path(root, "Code/sfa_functions.R"))

# read in data to run models on
db1 <- readRDS(file.path(root, "Cache/db1.rds"))
db1 <- unique(db1)
db1 <- filter(db1, lab < 2000, logN >= 0)

# make some new variables for convenience.
db1$logNloglab <- db1$logN * db1$loglab
db1$logslope <- log(db1$slope + 1)
db1$agesq <- db1$age^2
db1$elevationsqt <- sqrt(db1$elevation)
db1$SOC2sq <- db1$SOC2^2
db1$GGD <- db1$GGD/1000

# -------------------------------------
# translog production function with 
# endogeneity

# there are some missing values so we exclude these from the data now to avoid
# problems later
db1 <- as.data.frame(model.matrix(~ -1 + logyld +  logN + loglab  +
                                    logNloglab + Pn + dist_market + sex + age + agesq +
                                    ed_any + yesN + crop_count2 +
                                    logarea + impr + logslope + elevationsqt +
                                    SOC2 + SOC2sq + phdum55_2_70 + dumoxen + title +
                                    rain_wq + relprice + Pm + GGD +
                                    cost2large_town + extension, data=db1))

# basic frontier model for comparison
sf <- summary(sfa(logyld ~ logN + loglab + logNloglab + logarea +
                    logslope + elevationsqt + crop_count2 + SOC2 +
                    phdum55_2_70 + GGD, data=db1),
      extraPar = TRUE)

# basic reduced model for comparison
RM <- lm(logN ~ relprice + age + agesq + dist_market +
           logslope + elevationsqt + crop_count2 + SOC2 +
           phdum55_2_70 + GGD, data=db1)

# data
# outcome of structural equation
Y <- db1$logyld

# exogenous variables present in reduced and
# structural models
X1 <- model.matrix(~ -1 + logslope + elevationsqt + crop_count2 + SOC2 +
                     phdum55_2_70 + GGD, data=db1)

# endogenous variable(s)
X2 <- model.matrix( ~ -1 + logN, data=db1)

# variables we want in the structural equation but not 
# the reduced equation
X3 <- model.matrix(~ -1 + loglab + logNloglab + logarea, data = db1)

# full X ( all structural variables )
X <- cbind(X2, X3, X1)

# instruments/variables only in reduced equation
W <- model.matrix(~ -1 + relprice + age + agesq + dist_market,
                  data = db1)

# complete set of regressors for reduced equation
Z <- cbind(1, W, X1)

# parameters for optimisation
n <- (ncol(X) + ncol(Z) + 5)
pars <- rep(1, n)

# test out the function based on
# starting parameters - starting parameters
# must give a non NAN answer
liml1(pars, X, X2, Y, Z)

# optimx function
stats <- optimx(pars, liml1, NULL,
                hess=NULL,
                lower=-Inf,
                upper=Inf,
                method="BFGS",
                itnmax=NULL,
                hessian=FALSE,
                control = list(trace=2, maxit=200),
                X, X2, Y, Z)

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
c(alpha, beta) # parameters from liml second stage
dimnames(X)[[2]] # order in which parameters appear from liml output
sf$mleParam[, 1] # normal (no endogeneity) SF ML estimates
sf$olsParam[, 1] # ols estimates
c(sigma2v, sigma2u) # variance and relative variance estimates

# match up Pi parameters with the names from the Z
# matrix (first stage)
Pi # parameters from liml first stage
dimnames(Z)[[2]] # order in which parameters appear from liml output
RM # ols normal reduced model

