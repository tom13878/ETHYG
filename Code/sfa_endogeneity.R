#'==============================================================================
#' Project:  IMAGINE ETH
#' Subject:  sfa with endogeneity
#' Author:   Michiel van Dijk & Tom Morley
#' Contact:  michiel.vandijk@wur.nl, Tomas.morley@wur.nl
#' Output:   Results table of sfa endogeneity 
#'==============================================================================

# get packages
library(pacman)
p_load(char=c("rprojroot", "boot", "dplyr", "frontier"), install=TRUE)

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
                                    cost2large_town, data=db1))

# basic frontier model
sf1 <- summary(sfa(logyld ~ logN + loglab  + logarea + crop_count2 + 
      logslope + impr + SOC2 + GGD + phdum55_2_70 + dumoxen, data=db1),
      extraPar = TRUE)

# data
# outcome of structural equation
Y <- db1$logyld

# exogenous variables that we want in reduced
# and structural equation
X1 <- model.matrix( ~ -1 + crop_count2 + 
                      logslope + impr + SOC2 + GGD + phdum55_2_70, data = db1)

# endogenous variable(s)
X2 <- model.matrix( ~ -1 + logN, data=db1)

# variables we want in the structural equation but not 
# the reduced equation
X3 <- model.matrix(~ -1 + loglab + logNloglab + logarea, data = db1)

# full X ( all structural variables )
X <- cbind(X1, X2, X3)

# instruments/variables only in reduced equation
W <- model.matrix(~ -1 + relprice + sex + age + agesq + ed_any + dist_market + cost2large_town,
                  data = db1)

# complete set of regressors for reduced equation
Z <- cbind(1, W, X1)

# initial parameters
pars <- c(5, rep(0.2, ncol(X)), 0.5, 2, 0.1, 0.6, 4, rep(0.1, ncol(Z)-1))

# try out the function
liml1(pars, X, X2, Y, Z)

# optimize function
stats <- optim(pars, liml1, NULL, X, X2, Y, Z, method="BFGS", control = list(trace=2, maxit=200),hessian=T)

# get parameters (same order as they appear in the liml1 function)
alpha <- stats$par[1]
beta <- stats$par[2:(1 + ncol(X))]
sigma2v <- stats$par[(2 + ncol(X))]
sigma2u <- stats$par[(3 + ncol(X))]
pos <- (4+ncol(X))
Sigmavn <- stats$par[pos]
Sigmann <- stats$par[pos + 1]
Pi <- stats$par[(pos + 2): length(stats$par)]

# match up beta and Pi parameters with the names from
# the X and Z matrices
beta
dimnames(X)[[2]]

Pi
dimnames(Z)[[2]]

# get the SEs from the hessian
SEs <- sqrt(diag(solve(stats$hessian)))

