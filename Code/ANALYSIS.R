# use full model to calculate the MPP etc

library(pacman)
p_load(c("frontier", "rprojroot"))
root <- find_root(is_rstudio_project)

# get data
db1 <- readRDS(file.path(root, "Cache/db1.rds"))
db1 <- unique(db1)

# there is no predict function in the sfa package
# but we would like one - make our own and source
# it in
source(file.path(root, "Code/predictsfa.R"))

# get square and interaction terms for translog
db1$logNsq <- db1$logN^2
db1$loglabsq <- db1$loglab^2
db1$logseedsq <- db1$logseed^2
db1$logarea_tot <- log(db1$area_tot)

# Also divide GGD and AI by 1000
# to get a more interpretable coef
db1$GGD <- db1$GGD/1000
db1$AI <- db1$AI/1000

#
sf11x9 <- sfa(logyld ~ logN + loglab + logseed +
                logNsq + loglabsq + logseedsq +
                logN:loglab + logN:logseed +
                loglab:logseed + logarea + phdum55_2_70 +
                crop_count2 + dumoxen + SOC2 + logslope +
                elevation + GGD + AI + TS|
                -1 + age + sex + ed_any + title +
                extension + credit + dist_market +
                popEA + logarea_tot, data=db1)

# we want to evaluate the frontier function
# and avoid the Z variables
zidx <- grep("Z_", names(coef(sf11x9)))
xvars <- names(coef(sf11x9))[1:(min(zidx)-1)]
X <- sf11x9$dataTable[, xvars]
X <- as.data.frame(X)
xcoef <- coef(sf11x9)[1:(min(zidx)-1)]
relprices <- db1$relprice[sf11x9$validObs]

# function to calculate the MPP
calc_mpp <- function(N, row){
  row$logN <- log(N)
  logY <- as.matrix(row) %*% xcoef
  Y <- exp(logY)
  MPP <- with(row, ((xcoef["logN"] + 
                       2*xcoef["logNsq"]*logN +
                       xcoef["logN:loglab"]*loglab +
                       xcoef["logN:logseed"]*logseed )*
                      (Y/N)))
  MPP
}

# function to calculate the mpp less
# the relative price 
mpp_price <- function(N, row, relprice){
  calc_mpp(N, row) - relprice
}

# function to calculate the economically
# optimal level of nitrogen
econ_opt <- function(i){
  
  # get a plot observation and relative
  # price
  row <- X[i, ]
  relprice <- relprices[i]
  
  # first we need a point above zero
  # but before the root to act as the
  # lower limit of the root finding 
  # function - We cap at 700 because
  # even on the most productive farms
  # 700 kg/ha of nitrogen is excessive
  lower <- tryCatch(optimize(function(x) mpp_price(x, row, relprice),
                       interval=c(1, 700),
                       maximum =TRUE)$maximum, error=function(err) NA)
  
  # then we need to find the root
  # i.e. the point at which the
  # nitrogen level crosses zero. 
  # this is the economically
  # constrained level of nitrogen
  tryCatch(uniroot(function(x) {mpp_price(x, row, relprice)},
                   interval=c(lower, 700))$root,
           error=function(err) NA)
}

# calculate the economically
# optimal level of nitrogen per plot
db2 <- db1[sf11x9$validObs, ]
db2$Npm <- sapply(1:nrow(X), econ_opt)

# we can calculate a mpp wrt nitrogen
# at the observed nitrogen values.
# For values of N that are less than 1
# the MPP will be negative and large
# because the formula for the mpp of a
# translog involves the term 1/N
# for a value of N = 0.01, for example,
# this will result in an multiplicative
# of 100.
# of nitrogen that are 0, this will be 
# negative
db2$mpp <- calc_mpp(exp(X$logN), X)
db2$mpp[db2$mpp < 0] <- NA


model <- sf11x9
# 1. Technical efficiency yield is found using the
# output of the sfa model
# Observations where Npm cannot be calculated are removed
db3 <- db2 %>%
  rename(Y = yld) %>%
  mutate(
    Ycor = exp(as.numeric(fitted(model))+log(as.numeric(efficiencies(model)))), 
    err = Ycor-Y,
    TEY = exp(as.numeric(fitted(model))),
    TE = as.numeric(efficiencies(model)),
    resid = as.numeric(resid(model))
)

# 2. Economic yield is found by evaluating the frontier 
# function at the economically optimal nitrogen rate (Npm)
# but also need this for the 
X$logN <- log(db3$Npm) # swap N for econ optimal N
db3$EY <- exp(as.matrix(X) %*% xcoef)
X$logN <- db3$logN # swap back

# 3. PFY: Feasible yield
# evaluate frontier function at N = 150
# Increase labour and seed rate by 10%
# turn on all dummies
X$logN <- log(150) # swap N for econ optimal N
X$log
db3$EY <- exp(as.matrix(X) %*% xcoef)
X$logN <- db3$logN # swap back

group_by(db2, ZONE) %>%
  summarise(n = n(),
            Nopt = mean(Npm, na.rm=TRUE),
            n2 = sum(N > 0),
            MPP = mean(mpp, na.rm=TRUE))
