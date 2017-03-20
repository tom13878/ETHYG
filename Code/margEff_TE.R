#'==============================================================================
#' Project:  IMAGINE ETH
#' Subject:  Exogenous marginal effects calcualtion
#' Author:   Michiel van Dijk & Tom Morley
#' Contact:  michiel.vandijk@wur.nl, tomas.morley@wur.nl
#' Output:   Blueprint for marginal effects of exogenous variables
#'==============================================================================

# get packages
library(pacman)
p_load(char=c("rprojroot", "optimx", "frontier", "truncnorm", "abind", "dplyr"), install=TRUE)

# get root path
root <- find_root(is_rstudio_project)

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

# get ML parameters using the frontier package
sf <- sfa(logyld ~ logN + loglab | -1 +
            extension + title + age,
          data=db1)
sfsum <- summary(sf, extraPar = TRUE)$mleParam[, 1]

# Get ML parameters from optimization output
xvars <- c("logN", "loglab")
zvars <- c("extension", "title", "age")
zcoef <- grep("Z_", names(sfsum))
beta <- sfsum[1:(min(zcoef)-1)]
sigma2u <- sfsum["sigmaSqU"]
sigma2v <- sfsum["sigmaSqV"]
delta <- sfsum[zcoef]

# get model matrix for Z coefficients 
Z <- sf$dataTable[, zvars]
Y <- sf$dataTable[, "logyld"]
X <- sf$dataTable[, xvars]

# Derive values from ML estimates (this will also work
# with frontier package sfa function output)
sigma2 <- sigma2u + sigma2v
sigma <- sqrt(sigma2)
mu <- Z %*% delta
epsilon <- residuals(sf)
mus <- (mu * sigma2v - epsilon * sigma2u)/sigma2
sigma2s <- sigma2u*sigma2v/sigma2
sigmas <- sqrt(sigma2s)
m <- mus/sigmas
g <- dnorm(m)/pnorm(m)

# calcualte marginal effects following method
# of Kumbhakar and Sun (2013)
scale <- as.numeric((sigma2v/sigma2) * (1 - m * g - g^2))
margEff <- outer(scale, delta)

# calculate the APE
APE <- apply(margEff, 2, mean)

# Kumbhakar and Sun (2015) recommend bootstrapping
# the procedure in order to build a confidence
# interval around the marginal effects estimates
# however, we should be aware that the marginal
# effects are observation specific. This is a 
# parametric bootstrap -> no resampling the data
# instead we draw from the distributions of u and
# v (could make our own truncated normal sampler
# using accept-reject sampling)

vi <- rnorm(length(Y), mean=0, sqrt(sigma2v))
ui <- rtruncnorm(length(Y), a = 0, b = Inf, mean=mu, sqrt(sigma2u))

# generate a new Y values based on the new vi and ui terms
logyld <- cbind(1, X) %*% beta + vi - ui

# using {Ynew, X, Z} as before, get estimated
# parameters and then calculate marginal effects
new_dat <- data.frame(logyld, X, Z)

# get ML parameters using the frontier package
sf <- sfa(logyld ~ logN + loglab | -1 +
            extension + title + age,
          data=new_dat)
sfsum <- summary(sf, extraPar = TRUE)$mleParam[, 1]

# Get ML parameters from optimization output
beta2 <- sfsum[1:(min(zcoef)-1)]
sigma2u2 <- sfsum["sigmaSqU"]
sigma2v2 <- sfsum["sigmaSqV"]
delta2 <- sfsum[zcoef]

# Derive values from ML estimates (this will also work
# with frontier package sfa function output)
sigma22 <- sigma2u2 + sigma2v2
sigma2 <- sqrt(sigma22)
mu2 <- Z %*% delta2
epsilon2 <- residuals(sf)
mus2 <- (mu2 * sigma2v2 - epsilon2 * sigma2u2)/sigma22
sigma2s2 <- sigma2u2*sigma2v2/sigma22
sigmas2 <- sqrt(sigma2s2)
m2 <- mus2/sigmas2
g2 <- dnorm(m2)/pnorm(m2)

# calcualte marginal effects following method
# of Kumbhakar and Sun (2013)
scale2 <- as.numeric((sigma2v2/sigma22) * (1 - m2 * g2 - g2^2))
margEff2 <- outer(scale2, delta2)

# and then basically repeat this procedure B
# times. Of course we cannot hold this information
# in a 2 by 2 matrix because we have N observations
# for p exogenous variables over B bootstrap samples
# so instead we use a three dimensional array

B <- 10
N <- length(Y)
p <- ncol(Z)
margEff_array <- array(data=NA, dim=c(N, p, B))

system.time({
# bootstrap
for (i in 1:B){
  vi <- rnorm(length(Y), mean=0, sqrt(sigma2v))
  ui <- rtruncnorm(length(Y), a = 0, b = Inf, mean=mu, sqrt(sigma2u))
  
  # generate a new Y values based on the new vi and ui terms
  logyld <- cbind(1, X) %*% beta + vi - ui
  
  # using {Ynew, X, Z} as before, get estimated
  # parameters and then calculate marginal effects
  new_dat <- data.frame(logyld, X, Z)
  
  # get ML parameters using the frontier package
  sf <- sfa(logyld ~ logN + loglab | -1 +
              extension + title + age,
            data=new_dat)
  sfsum <- summary(sf, extraPar = TRUE)$mleParam[, 1]
  
  # Get ML parameters from optimization output
  beta2 <- sfsum[1:(min(zcoef)-1)]
  sigma2u2 <- sfsum["sigmaSqU"]
  sigma2v2 <- sfsum["sigmaSqV"]
  delta2 <- sfsum[zcoef]
  
  # Derive values from ML estimates (this will also work
  # with frontier package sfa function output)
  sigma22 <- sigma2u2 + sigma2v2
  sigma2 <- sqrt(sigma22)
  mu2 <- Z %*% delta2
  epsilon2 <- residuals(sf)
  mus2 <- (mu2 * sigma2v2 - epsilon2 * sigma2u2)/sigma22
  sigma2s2 <- sigma2u2*sigma2v2/sigma22
  sigmas2 <- sqrt(sigma2s2)
  m2 <- mus2/sigmas2
  g2 <- dnorm(m2)/pnorm(m2)
  
  # calcualte marginal effects following method
  # of Kumbhakar and Sun (2013)
  scale2 <- as.numeric((sigma2v2/sigma22) * (1 - m2 * g2 - g2^2))
  margEff_array[, c(1:ncol(Z)), i] <- outer(scale2, delta2)
}
})
# bind the original margEffs to the parametric bootstrap
# margEffs in order to get sample of (1 + B) marginal
# effects for each observation
margEff_array <- abind(margEff_array, margEff, along=3)

# now for each observation we have B + 1 parametric
# bootstrap samples (including the first one)
# for each variable. From this we can calculate the
# 95% confidence interval.

# for example for the first observation, for the
# variable extension the 95% confidence interval
# would be

quantile(margEff_array[1,1,], probs=c(0.025, 0.975))

# and the actual estimate of the marginal effects would
# be
margEff[1, 1]

# we can do this for all observations/variables
# so that the confidence intervals are now stored
# in the third dimension
CI <- apply(margEff_array, c(1, 2), function(x) quantile(x, probs=c(0.025, 0.975)))
CI[, 1, 1] # check against quantile function above

# given we have an APE for each observation
# we could do something similar for the CI
APE_CI <- apply(CI, c(1, 3), mean)

# check APE falls inside CI
APE_CI
APE 
