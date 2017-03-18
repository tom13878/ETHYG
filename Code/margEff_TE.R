#'==============================================================================
#' Project:  IMAGINE ETH
#' Subject:  Exogenous marginal effects calcualtion
#' Author:   Michiel van Dijk & Tom Morley
#' Contact:  michiel.vandijk@wur.nl, Tomas.morley@wur.nl
#' Output:   Blueprint for marginal effects of exogenous variables
#'==============================================================================

# get packages
library(pacman)
p_load(char=c("optimx", "frontier", "truncnorm", "abind"), install=TRUE)

# include exogenous determinants of inefficiency
# basic SF log likelihood function
SFlogLik_x <- function(pars, X, Y, Z){
  
  # get parameters
  alpha <- pars[1]
  beta <- pars[2:(1 + ncol(X))]
  sigma2u <- pars[(2 + ncol(X))]
  sigma2v <- pars[(3 + ncol(X))]
  delta <- pars[(4 + ncol(X)):length(pars)]
  
  # create variables
  epsilon <- Y - alpha - (X %*% beta)
  mu <- Z %*% delta
  sigma2 <- sigma2u + sigma2v
  sigma <- sqrt(sigma2)
  sigmau <- sqrt(sigma2u)
  sigmav <- sqrt(sigma2v)
  lambda <- sigmau/sigmav
  mus <- (mu * sigma2v - epsilon * sigma2u)/sigma2
  sigma2s <- sigma2u*sigma2v/sigma2
  sigmas <- sqrt(sigma2s)
  
  # construct (log) densities
  Li <- -(1/2) * log(sigma2) + dnorm((mu + epsilon)/sigma, log=TRUE) +
    pnorm(mus/sigmas, log=TRUE) - pnorm(mu/sigmau, log=TRUE)
  
  
  # evaluate log density and sum
  -sum(Li)
}

dat <- model.matrix(~ -1 + logyld + logN + logseed + loglab +
                      extension + title, data=panel_ETH)
yindx <- which(dimnames(dat)[[2]] == "logyld")
zindx <- which(dimnames(dat)[[2]] == c("extension", "title"))

Y <- dat[, yindx]
Z <- dat[, zindx]
X <- dat[, -c(yindx, zindx)]

# initial parameters
pars <- rep(0.5, ncol(X) + ncol(Z) + 3)

# test function works at initial parameters
# if NA is returned, pick better starting
# parameters
SFlogLik_x(pars, X, Y, Z)

# use optimx to get ML estimates
stats <- optimx(pars, SFlogLik_x, NULL,
                hess=NULL,
                lower=-Inf,
                upper=Inf,
                method="BFGS",
                itnmax=NULL,
                hessian=FALSE,
                control = list(trace=2, maxit=200),
                X, Y, Z)

# check if it matches output from frontier package
sf <- sfa(logyld ~ logN + logseed + loglab |
            -1 + extension + title, data=panel_ETH)

# Get ML parameters from optimization output
alpha <- as.numeric(stats[1])
beta <- as.numeric(stats[2:(1 + ncol(X))])
sigma2u <- as.numeric(stats[(2 + ncol(X))])
sigma2v <- as.numeric(stats[(3 + ncol(X))])
delta <- as.numeric(stats[(4 + ncol(X)):length(pars)])

# Derive values from ML estimates (this will also work
# with frontier package sfa function output)
sigma2 <- sigma2u + sigma2v
sigma <- sqrt(sigma2)
mu <- Z %*% delta
epsilon <- Y - alpha - (X %*% beta)
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
Ynew <- alpha + X %*% beta + vi - ui

# using {Ynew, X, Z} as before, get estimated
# parameters and then calculate marginal effects
get_stats <- function(Ynew){
  optimx(pars, SFlogLik_x, NULL,
                hess=NULL,
                lower=-Inf,
                upper=Inf,
                method="BFGS",
                itnmax=NULL,
                hessian=FALSE,
                control = list(trace=2, maxit=200),
                X, Ynew, Z)
}

# Get ML parameters from optimization output
alpha2 <- as.numeric(stats[1])
beta2 <- as.numeric(stats[2:(1 + ncol(X))])
sigma2u2 <- as.numeric(stats[(2 + ncol(X))])
sigma2v2 <- as.numeric(stats[(3 + ncol(X))])
delta2 <- as.numeric(stats[(4 + ncol(X)):length(pars)])

# Derive values from ML estimates (this will also work
# with frontier package sfa function output)
sigma22 <- sigma2u2 + sigma2v2
sigma2 <- sqrt(sigma22)
mu2 <- Z %*% delta2
epsilon2 <- Ynew - alpha2 - (X %*% beta2)
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

B <- 10
N <- length(Y)
p <- ncol(Z)
margEff_array <- array(data=NA, dim=c(N, p, B))

# bootstrap
for (i in 1:B){
  vi <- rnorm(length(Y), mean=0, sqrt(sigma2v))
  ui <- rtruncnorm(length(Y), a = 0, b = Inf, mean=mu, sqrt(sigma2u))
  Ynew <- alpha + X %*% beta + vi - ui 
  stats_new <- get_stats(Ynew)
  
  # Get ML parameters from optimization output
  alpha2 <- as.numeric(stats_new[1])
  beta2 <- as.numeric(stats_new[2:(1 + ncol(X))])
  sigma2u2 <- as.numeric(stats_new[(2 + ncol(X))])
  sigma2v2 <- as.numeric(stats_new[(3 + ncol(X))])
  delta2 <- as.numeric(stats_new[(4 + ncol(X)):length(pars)])
  
  # Derive values from ML estimates (this will also work
  # with frontier package sfa function output)
  sigma22 <- sigma2u2 + sigma2v2
  sigma2 <- sqrt(sigma22)
  mu2 <- Z %*% delta2
  epsilon2 <- Ynew - alpha2 - (X %*% beta2)
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
