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

dat <- model.matrix(~ -1 + logyld + logN + logseed + loglab, data=panel_ETH)
indx <- which(dimnames(dat)[[2]] == "logyld")
Y <- dat[, indx]
X <- dat[, -indx]

sfa(logyld ~ logN + logseed + loglab, data=panel_ETH)

pars <- rep(0.5, ncol(X) + 3)
SFlogLik(pars, X, Y)
stats <- optimx(pars, SFlogLik, NULL,
                hess=NULL,
                lower=-Inf,
                upper=Inf,
                method="BFGS",
                itnmax=NULL,
                hessian=FALSE,
                control = list(trace=2, maxit=200),
                X, Y)
stats

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

# check if it matches output from sfa package
sf <- sfa(logyld ~ logN + logseed + loglab | -1 + extension + title, data=panel_ETH)

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

# Get ML parameters from optimization output
alpha <- as.numeric(stats[1])
beta <- as.numeric(stats[2:(1 + ncol(X))])
sigma2u <- as.numeric(stats[(2 + ncol(X))])
sigma2v <- as.numeric(stats[(3 + ncol(X))])
delta <- as.numeric(stats[(4 + ncol(X)):length(pars)])

# Derive values from ML estimates
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
# effects are observation specific. As a result,
# so are the confidence intervals. An alternative
# way to bootstrap a distribution for the APE, may 
# be to simplt repeat the above calcualtion B
# times and record the APE each time.
