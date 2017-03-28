#'==============================================================================
#' Project:  IMAGINE ETH
#' Subject:  Analysis file
#' Author:   Michiel van Dijk & Tom Morley
#' Contact:  michiel.vandijk@wur.nl, tomas.morley@wur.nl
#' Output:   Full analysis file
#'==============================================================================

# set project root
library(dplyr)
library(rprojroot)
library(moments)
root <- find_root(is_rstudio_project)

# load packages
library(pacman)
p_load(frontier)

# get data
db1 <- readRDS(file.path(root, "Cache/db1.rds"))

# make a summary tab of variables for paper
sum_dat <- select(db1, yld, N, area, lab, seedha,
                   phdum55_2_70, crop_count2,
                  dumoxen, SOC2, logslope,
                  elevation, GGD, AI, TS)

Mean <- colMeans(sum_dat, na.rm=TRUE)
Median <- apply(sum_dat, 2, function(col) median(col, na.rm=TRUE))
SD <- apply(sum_dat, 2, function(col) sd(col, na.rm=TRUE))
Skewness <- apply(sum_dat, 2, function(col) skewness(col, na.rm=TRUE))
Min <- apply(sum_dat, 2, function(col) min(col, na.rm=TRUE))
Max <- apply(sum_dat, 2, function(col) max(col, na.rm=TRUE))
sum_tab <- data.frame(Variable=colnames(sum_dat),
                      Mean, Median, SD, Skewness,
                      Min, Max)
sum_tab[, -1] <- round(sum_tab[,-1], 3)
row.names(sum_tab) <- NULL

# basic translog function
sf1 <- sfa(logyld ~ logN + loglab + logseed +
            logNsq + loglabsq + logseedsq +
            logN:loglab + logN:logseed +
            loglab:logseed, data=db1)

# # translog function with logarea as a regressor
# sf2 <- sfa(logyld ~ logN + loglab + logseed +
#             logNsq + loglabsq + logseedsq +
#             logN:loglab + logN:logseed +
#             loglab:logseed + logarea
#             , data=db1)
# 
# # translog function with a dummy if the ph is
# # in the best range for maize
# sf3 <- sfa(logyld ~ logN + loglab + logseed +
#              logNsq + loglabsq + logseedsq +
#              logN:loglab + logN:logseed +
#              loglab:logseed + logarea + phdum55_2_70
#            , data=db1)
# 
# # translog function adding a dummy for
# # crop count
# sf4 <- sfa(logyld ~ logN + loglab + logseed +
#              logNsq + loglabsq + logseedsq +
#              logN:loglab + logN:logseed +
#              loglab:logseed + logarea + phdum55_2_70 +
#              crop_count2
#            , data=db1)
# 
# # translog function adding a dummy if the 
# # farmer owns oxen
# sf5 <- sfa(logyld ~ logN + loglab + logseed +
#              logNsq + loglabsq + logseedsq +
#              logN:loglab + logN:logseed +
#              loglab:logseed + logarea + phdum55_2_70 +
#              crop_count2 + dumoxen
#            , data=db1)
# 
# # translog function adding a variable for
# # soil organic content
# sf6 <- sfa(logyld ~ logN + loglab + logseed +
#              logNsq + loglabsq + logseedsq +
#              logN:loglab + logN:logseed +
#              loglab:logseed + logarea + phdum55_2_70 +
#              crop_count2 + dumoxen + SOC2 
#            , data=db1)
# 
# # translog function adding a variable for
# # the log of the plot slope. Note that
# # once we account for slope, the logarea 
# # variable becomes significant
# sf7 <- sfa(logyld ~ logN + loglab + logseed +
#              logNsq + loglabsq + logseedsq +
#              logN:loglab + logN:logseed +
#              loglab:logseed + logarea + phdum55_2_70 +
#              crop_count2 + dumoxen + SOC2 + logslope 
#            , data=db1)
# 
# # translog function adding a variable for
# # the elevation of the plot. Note the
# # effect on nitrogen
# sf8 <- sfa(logyld ~ logN + loglab + logseed +
#              logNsq + loglabsq + logseedsq +
#              logN:loglab + logN:logseed +
#              loglab:logseed + logarea + phdum55_2_70 +
#              crop_count2 + dumoxen + SOC2 + logslope +
#              elevation
#            , data=db1)
# 
# # translog function adding a variable for
# # the growing degree days. 
# sf9 <- sfa(logyld ~ logN + loglab + logseed +
#              logNsq + loglabsq + logseedsq +
#              logN:loglab + logN:logseed +
#              loglab:logseed + logarea + phdum55_2_70 +
#              crop_count2 + dumoxen + SOC2 + logslope +
#              elevation + GGD
#            , data=db1)
# 
# # translog function adding a variable for
# # the aridity index (AI). 
# sf10 <- sfa(logyld ~ logN + loglab + logseed +
#              logNsq + loglabsq + logseedsq +
#              logN:loglab + logN:logseed +
#              loglab:logseed + logarea + phdum55_2_70 +
#              crop_count2 + dumoxen + SOC2 + logslope +
#              elevation + GGD + AI
#            , data=db1)
# 
# # translog function adding a variable for
# # the temperate seasonality (TS). 
sf11 <- sfa(logyld ~ logN + loglab + logseed +
              logNsq + loglabsq + logseedsq +
              logN:loglab + logN:logseed +
              loglab:logseed + logarea + phdum55_2_70 +
              crop_count2 + dumoxen + SOC2 + logslope +
              elevation + GGD + AI + TS
            , data=db1)
# 
# summary(sf1)
# summary(sf11)

# now try with determinants of technical inefficiency
# using sf1
# sf1x1 <- sfa(logyld ~ logN + loglab + logseed +
#              logNsq + loglabsq + logseedsq +
#              logN:loglab + logN:logseed +
#              loglab:logseed |
#               -1 + age, data=db1)
# 
# sf1x2 <- sfa(logyld ~ logN + loglab + logseed +
#               logNsq + loglabsq + logseedsq +
#               logN:loglab + logN:logseed +
#               loglab:logseed |
#               -1 + age + sex, data=db1)
# 
# sf1x3 <- sfa(logyld ~ logN + loglab + logseed +
#                logNsq + loglabsq + logseedsq +
#                logN:loglab + logN:logseed +
#                loglab:logseed |
#                -1 + age + sex + ed_any, data=db1)
# 
# sf1x4 <- sfa(logyld ~ logN + loglab + logseed +
#                logNsq + loglabsq + logseedsq +
#                logN:loglab + logN:logseed +
#                loglab:logseed |
#                -1 + age + sex + ed_any + title, data=db1)
# 
# sf1x5 <- sfa(logyld ~ logN + loglab + logseed +
#                logNsq + loglabsq + logseedsq +
#                logN:loglab + logN:logseed +
#                loglab:logseed |
#                -1 + age + sex + ed_any + title +
#                extension, data=db1)
# 
# sf1x6 <- sfa(logyld ~ logN + loglab + logseed +
#                logNsq + loglabsq + logseedsq +
#                logN:loglab + logN:logseed +
#                loglab:logseed |
#                -1 + age + sex + ed_any + title +
#                extension + credit, data=db1)
# 
# sf1x7 <- sfa(logyld ~ logN + loglab + logseed +
#                logNsq + loglabsq + logseedsq +
#                logN:loglab + logN:logseed +
#                loglab:logseed |
#                -1 + age + sex + ed_any + title +
#                extension + credit + dist_market, data=db1)
# 
# sf1x8 <- sfa(logyld ~ logN + loglab + logseed +
#                logNsq + loglabsq + logseedsq +
#                logN:loglab + logN:logseed +
#                loglab:logseed |
#                -1 + age + sex + ed_any + title +
#                extension + credit + dist_market +
#                popEA, data=db1)

sf1x9 <- sfa(logyld ~ logN + loglab + logseed +
               logNsq + loglabsq + logseedsq +
               logN:loglab + logN:logseed +
               loglab:logseed |
               -1 + age + sex + ed_any + title +
               extension + credit + dist_market +
               popEA + logarea_tot, data=db1)

# now try with determinants of technical inefficiency
# using sf11
# sf11x1 <- sfa(logyld ~ logN + loglab + logseed +
#                logNsq + loglabsq + logseedsq +
#                logN:loglab + logN:logseed +
#                loglab:logseed + logarea + phdum55_2_70 +
#                 crop_count2 + dumoxen + SOC2 + logslope +
#                 elevation + GGD + AI + TS|
#                -1 + age, data=db1)
# 
# sf11x2 <- sfa(logyld ~ logN + loglab + logseed +
#                logNsq + loglabsq + logseedsq +
#                logN:loglab + logN:logseed +
#                loglab:logseed + logarea + phdum55_2_70 +
#                 crop_count2 + dumoxen + SOC2 + logslope +
#                 elevation + GGD + AI + TS|
#                -1 + age + sex, data=db1)
# 
# sf11x3 <- sfa(logyld ~ logN + loglab + logseed +
#                logNsq + loglabsq + logseedsq +
#                logN:loglab + logN:logseed +
#                loglab:logseed + logarea + phdum55_2_70 +
#                 crop_count2 + dumoxen + SOC2 + logslope +
#                 elevation + GGD + AI + TS|
#                -1 + age + sex + ed_any, data=db1)
# 
# sf11x4 <- sfa(logyld ~ logN + loglab + logseed +
#                logNsq + loglabsq + logseedsq +
#                logN:loglab + logN:logseed +
#                loglab:logseed + logarea + phdum55_2_70 +
#                 crop_count2 + dumoxen + SOC2 + logslope +
#                 elevation + GGD + AI + TS|
#                -1 + age + sex + ed_any + title, data=db1)
# 
# sf11x5 <- sfa(logyld ~ logN + loglab + logseed +
#                logNsq + loglabsq + logseedsq +
#                logN:loglab + logN:logseed +
#                loglab:logseed + logarea + phdum55_2_70 +
#                 crop_count2 + dumoxen + SOC2 + logslope +
#                 elevation + GGD + AI + TS|
#                -1 + age + sex + ed_any + title +
#                extension, data=db1)
# 
# sf11x6 <- sfa(logyld ~ logN + loglab + logseed +
#                logNsq + loglabsq + logseedsq +
#                logN:loglab + logN:logseed +
#                loglab:logseed + logarea + phdum55_2_70 +
#                 crop_count2 + dumoxen + SOC2 + logslope +
#                 elevation + GGD + AI + TS|
#                -1 + age + sex + ed_any + title +
#                extension + credit, data=db1)
# 
# sf11x7 <- sfa(logyld ~ logN + loglab + logseed +
#                logNsq + loglabsq + logseedsq +
#                logN:loglab + logN:logseed +
#                loglab:logseed + logarea + phdum55_2_70 +
#                 crop_count2 + dumoxen + SOC2 + logslope +
#                 elevation + GGD + AI + TS|
#                -1 + age + sex + ed_any + title +
#                extension + credit + dist_market, data=db1)
# 
# sf11x8 <- sfa(logyld ~ logN + loglab + logseed +
#                logNsq + loglabsq + logseedsq +
#                logN:loglab + logN:logseed +
#                loglab:logseed + logarea + phdum55_2_70 +
#                 crop_count2 + dumoxen + SOC2 + logslope +
#                 elevation + GGD + AI + TS|
#                -1 + age + sex + ed_any + title +
#                extension + credit + dist_market +
#                popEA, data=db1)

sf11x9 <- sfa(logyld ~ logN + loglab + logseed +
               logNsq + loglabsq + logseedsq +
               logN:loglab + logN:logseed +
               loglab:logseed + logarea + phdum55_2_70 +
                crop_count2 + dumoxen + SOC2 + logslope +
                elevation + GGD + AI + TS|
               -1 + age + sex + ed_any + title +
               extension + credit + dist_market +
               popEA + logarea_tot, data=db1)

# make a table of the results from
# the translog production and the
# translog production with the 
# exogenous variables

model1 <- sprintf("%.3f", round(summary(sf1)$mleParam[, 1], 3))
model2 <- sprintf("%.3f", round(summary(sf11)$mleParam[, 1], 3))
model3 <- sprintf("%.3f", round(summary(sf1x9)$mleParam[, 1], 3))
model4 <- sprintf("%.3f", round(summary(sf11x9)$mleParam[, 1], 3))

# get stars from pvals
p1 <- cut(summary(sf1)$mleParam[, 4], breaks=c(0, 0.01, 0.05, 1), labels=c("**", "*", ""), include.lowest = TRUE)
p2 <- cut(summary(sf11)$mleParam[, 4], breaks=c(0, 0.01, 0.05, 1), labels=c("**", "*", ""), include.lowest = TRUE)
p3 <- cut(summary(sf1x9)$mleParam[, 4], breaks=c(0, 0.01, 0.05, 1), labels=c("**", "*", ""), include.lowest = TRUE)
p4 <- cut(summary(sf11x9)$mleParam[, 4], breaks=c(0, 0.01, 0.05, 1), labels=c("**", "*", ""), include.lowest = TRUE)

# add stars to models
model1 <- as.data.frame(paste(model1, p1, sep=""))
model2 <- as.data.frame(paste(model2, p2, sep=""))
model3 <- as.data.frame(paste(model3, p3, sep=""))
model4 <- as.data.frame(paste(model4, p4, sep=""))

# set row names
model1$Variable <- names(coef(sf1))
model2$Variable <- names(coef(sf11))
model3$Variable <- names(coef(sf1x9))
model4$Variable <- names(coef(sf11x9))

# join all results
results_tab <- full_join(model1, model2) %>%
  full_join(model3) %>%
  full_join(model4) %>%
  select(Variable, everything())
names(results_tab) <- c("Variable", "model 1", "model 2",
                        "model 3", "model 4")
move <- which(results_tab$Variable %in% c("sigmaSq", "gamma"))
results_tab <- rbind(results_tab[-move, ], results_tab[move, ])
row.names(results_tab) <- NULL

# from the models which have exogenous determinants
# of inefficiency, calculate the marginal effects
# of each Z variable on the inefficiency

# get the location of the z variables
zidx <- grep("Z_", names(coef(sf1x9)))
zvars <- gsub("Z_", "", names(coef(sf1x9))[zidx])

# get ML parameters from output
beta <- summary(sf1x9)$mleParam[1:(min(zidx)-1)]
sigma2u <- summary(sf1x9, extraPar=TRUE)$mleParam[, 1]["sigmaSqU"]
sigma2v <- summary(sf1x9, extraPar=TRUE)$mleParam[, 1]["sigmaSqV"]
delta <- summary(sf1x9)$mleParam[zidx]

# get the Z matrix
Z <- sf1x9$dataTable[, zvars]

# Derive values from ML estimates 
sigma2 <- sigma2u + sigma2v
sigma <- sqrt(sigma2)
mu <- Z %*% delta
epsilon <- residuals(sf1x9)
mus <- (mu * sigma2v - epsilon * sigma2u)/sigma2
sigma2s <- sigma2u*sigma2v/sigma2
sigmas <- sqrt(sigma2s)
m <- mus/sigmas
g <- dnorm(m)/pnorm(m)

# calcualte marginal effects following method
# of Kumbhakar and Sun (2013)
scale <- as.numeric((sigma2v/sigma2) * (1 - m * g - g^2))
margEff <- outer(scale, delta)

# calculate the APE, to give a 
# consise idea of the marginal
# elasticity of inefficiency
# with regard to the z vars
APE_sf1x9 <- apply(margEff, 2, mean)

# do the same for the translog with all
# of the control variables.

# get the location of the z variables
zidx <- grep("Z_", names(coef(sf11x9)))
zvars <- gsub("Z_", "", names(coef(sf11x9))[zidx])

# get ML parameters from output
beta <- summary(sf11x9)$mleParam[1:(min(zidx)-1)]
sigma2u <- summary(sf11x9, extraPar=TRUE)$mleParam[, 1]["sigmaSqU"]
sigma2v <- summary(sf11x9, extraPar=TRUE)$mleParam[, 1]["sigmaSqV"]
delta <- summary(sf11x9)$mleParam[zidx]

# get the Z matrix
Z <- sf11x9$dataTable[, zvars]

# Derive values from ML estimates 
sigma2 <- sigma2u + sigma2v
sigma <- sqrt(sigma2)
mu <- Z %*% delta
epsilon <- residuals(sf11x9)
mus <- (mu * sigma2v - epsilon * sigma2u)/sigma2
sigma2s <- sigma2u*sigma2v/sigma2
sigmas <- sqrt(sigma2s)
m <- mus/sigmas
g <- dnorm(m)/pnorm(m)

# calcualte marginal effects following method
# of Kumbhakar and Sun (2013)
scale <- as.numeric((sigma2v/sigma2) * (1 - m * g - g^2))
margEff <- outer(scale, delta)

# calculate the APE, to give a 
# consise idea of the marginal
# elasticity of inefficiency
# with regard to the z vars
APE_sf11x9 <- apply(margEff, 2, mean)

# join APEs together to make a table
ME_tab <- round(data.frame(cbind(APE_sf1x9, APE_sf11x9)), 5)
ME_tab <- as.data.frame(lapply(ME_tab, function(col) sprintf("%.3f", round(col, 3))))
row.names(ME_tab) <- zvars
names(ME_tab) <- c("model 3", "model 4")

# take out trash
rm(APE_sf11x9, APE_sf1x9, beta, db1, delta, epsilon,
   g, m, margEff, mu, mus, scale,
   sigma, sigma2, sigma2s, sigma2u, sigma2v,
   sigmas, Z, zidx, zvars, sf1, sf11, sf11x9,
   sf1x9, model1, model2, model3, model4)
