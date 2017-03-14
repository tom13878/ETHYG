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

# basic frontier model for comparison
sf <- summary(sfa(logyld ~ logN + loglab + logNloglab + logarea +
                    logslope + elevationsqt + crop_count2 + SOC2 +
                    phdum55_2_70 + GGD, data=db1),
      extraPar = TRUE)

# basic reduced model for comparison
RM <- lm(logN ~ relprice + age + agesq + dist_market +
           logslope + elevationsqt + crop_count2 + SOC2 +
           phdum55_2_70 + GGD, data=db1)

# starting value for Sigmann
Sigmann_hat <- mean(residuals(RM)^2)

# simple ols model with translog for comparison
TL <- lm(logyld ~ logN + loglab + logNloglab + logarea +
           logslope + elevationsqt + crop_count2 + SOC2 +
           phdum55_2_70 + GGD, data=db1)

# starting value for Sigmavn
Sigmavn_hat <- mean(residuals(RM)*residuals(TL))

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
#  + sex + age + agesq + ed_any + dist_market + cost2large_town
W <- model.matrix(~ -1 + relprice + age + agesq + dist_market,
                  data = db1)

# complete set of regressors for reduced equation
Z <- cbind(1, W, X1)

# initial parameters - starting values were chosen
# on the basis of the output from simple models
# bearing in mind that with many parameters the 
# chances of falling into a local minima are high ->
# therefore, we would like to start our search
# for the correct parameters in the region of 
# parameter space that we would expect to find
# them
pars <- c(6, rep(0.2, ncol(X)), 0.4, 1.7, 0.03, 2.7, 7.5, rep(0.1, ncol(Z)-1))
# #pars[ncol(X) + 2] <- sfco["sigmaSqV"]
# #pars[ncol(X) + 3] <- sfco["sigmaSqU"]
# pars[ncol(X) + 4] <- Sigmavn_hat
# pars[ncol(X) + 5] <- Sigmann_hat
# try out the function
liml1(pars, X, X2, Y, Z)

# optimize function
stats <- optim(pars, liml1, NULL, X, X2, Y, Z,
               method="BFGS",
               control = list(trace=2, maxit=200),
               hessian=T)

# get parameters (same order as they appear in the liml1 function)
alpha <- stats$par[1]
beta <- stats$par[2:(1 + ncol(X))]
sigma2v <- stats$par[(2 + ncol(X))]
sigma2u <- stats$par[(3 + ncol(X))]
lambda <- sqrt(exp(sigma2u))/sqrt(exp(sigma2v))

pos <- (4+ncol(X))
Sigmavn <- stats$par[pos]
Sigmann <- stats$par[pos + 1]
Pi <- stats$par[(pos + 2): length(stats$par)]

# match up beta and Pi parameters with the names from
# the X and Z matrices
alpha
beta
dimnames(X)[[2]]
sf

Pi
dimnames(Z)[[2]]
RM

