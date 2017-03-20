#'==============================================================================
#' Project:  IMAGINE ETH
#' Subject:  sfa functions
#' Author:   Michiel van Dijk & Tom Morley
#' Contact:  michiel.vandijk@wur.nl, Tomas.morley@wur.nl
#' Output:   sfa functions 
#'==============================================================================

# basic SF log likelihood function
SFlogLik <- function(pars, X, Y){
  
  # get parameters
  alpha <- pars[1]
  beta <- pars[2:(1 + ncol(X))]
  sigma2u <- pars[(2 + ncol(X))]
  sigma2v <- pars[(3 + ncol(X))]
  
  # create variables
  epsilon <- Y - alpha - (X %*% beta)
  sigma2 <- sigma2u + sigma2v
  sigma <- sqrt(sigma2)
  sigmau <- sqrt(sigma2u)
  sigmav <- sqrt(sigma2v)
  lambda <- sigmau/sigmav
  
  # construct (log) densities
  f_epsilon <- log(2) - log(sigma) + dnorm(epsilon/sigma, log=TRUE) + pnorm(-lambda * epsilon/sigma, log=TRUE)
  
  # evaluate log density and sum
  -sum(f_epsilon)
}

# LIML 2 using the derivation in APS 2016 
liml2 <- function(pars, X, X2, Y, Z, n){
  
  # get all parameters
  alpha <- pars[1]
  beta <- pars[2:(1 + ncol(X))]
  sigma2v <- pars[(2 + ncol(X))]
  sigma2u <- pars[(3 + ncol(X))]
  pos <- (4+ncol(X))
  Sigmavn <- pars[pos]
  Sigmann <- pars[pos + 1]
  Pi <- pars[(pos + 2): length(pars)]
  
  # calculate values derived from parameters
  epsilon <- Y - alpha - (X %*% beta)
  eta <- X2 - (Z %*% Pi)
  muc <- Sigmavn * (1/Sigmann) * eta
  sigmau <- sqrt(sigma2u)
  sigma2c <- sigma2v - Sigmavn * (1/Sigmann) * Sigmavn
  sigmac <- sqrt(sigma2c)
  lambda <- sigmau/sigmac
  sigma2 <- sigma2u + sigma2c
  sigma <- sqrt(sigma2)
  
  # get first part of log likelihood
  lnL1 <- -(n/2) * log(sigma2) - (1/(2*sigma2)) * sum((epsilon - muc)^2) +
    sum(pnorm(-lambda * (epsilon - muc)/sigma, log = TRUE))
  
  # get second part of log likelihood
  lnL2 <- -(n/2) * log(Sigmann) - (1/2) * sum(eta * (1/Sigmann) * eta)
  
  -(lnL1 + lnL2)
  
}



