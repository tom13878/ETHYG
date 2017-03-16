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


# LIML for endogeneity
liml1 <- function(pars, X, X2, Y, Z){
  
  # get all parameters
  alpha <- pars[1]
  beta <- pars[2:(1 + ncol(X))]
  sigma2v <- exp(pars[(2 + ncol(X))]) + 10^-10
  sigma2u <- exp(pars[(3 + ncol(X))]) + 10^-10
  pos <- (4+ncol(X))
  Sigmavn <- pars[pos]
  Sigmann <- exp(pars[pos + 1]) + 10^-10
  Pi <- pars[(pos + 2): length(pars)]
  
  # calculate values derived from these
  epsilon <- Y - alpha - (X %*% beta)
  eta <- X2 - (Z %*% Pi)
  muc <- Sigmavn * (1/Sigmann) * eta
  sigmau <- sqrt(sigma2u)
  sigma2c <- exp(sigma2v - Sigmavn * (1/Sigmann) * Sigmavn) + 10^-10
  sigmac <- sqrt(sigma2c)
  lambda <- sigmau/sigmac
  sigma2 <- sigma2u + sigma2c
  sigma <- sqrt(sigma2)
  
  # construct (log) densities
  f_epsilon_eta <- -(1/2) * log(Sigmann) - (1/2) * eta * (1/Sigmann) * eta - log(sigma) +
    dnorm((epsilon - muc)/sigma, log=TRUE) + pnorm(-lambda * (eta - muc)/sigma, log=TRUE)
  
  # return the (minus) the sum of the log densities
  -sum(f_epsilon_eta)
}
