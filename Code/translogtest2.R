# set project root
library(rprojroot)
root <- find_root(is_rstudio_project)

# load packages
library(frontier)
library(moments)
library(dplyr)

# source in prepared data and functions
source(file.path(root, "Code/ETH_2013_prepare4analysis.R"))
source(file.path(root, "Code/translog_formula.R"))

# define inputs for translog
# define inputs for environmental variables
# define z inputs for regressing against technical
# efficiency

# Cobb Douglas ols
olsCD <- lm(logyld ~ noN + logN + loglab + 
              logarea +
              impr +
              slope + elevation +
              rain_wq + 
              crop_count2,
            data = db1)
summary(olsCD)

# Cobb Douglas sfa
sfaCD <- sfa(logyld ~ noN + logN + loglab + 
              logarea +
              impr +
              slope + elevation +
              rain_wq + 
              crop_count2,
            data = db1)
summary(sfaCD, extraPar = TRUE)
lrtest(sfaCD)

# Compare CD ols and sfa
# olsCD_CRE
model <- olsCD

# Note that MPP cannot be calculated for plots with N=0 and are therefore set to 0.
db_olsCD <- db1 %>% mutate(elastfert  = coef(model)["logN"],
                               imprdf = (coef(model)["impr"]*impr),
                               lnA = coef(model)["(Intercept)"] +
                                 (coef(model)["noN"]*noN) +
                                 (coef(model)["logarea"]*logarea) +
                                 (coef(model)["slope"]*slope) +
                                 (coef(model)["elevation"]*elevation) +
                                 (coef(model)["rain_wq"]* rain_wq) +
                                 (coef(model)["crop_count2"]*crop_count2),
                               lnA2 = lnA + imprdf,
                               constantfactor = exp(lnA2)*elastfert*(lab^coef(model)["loglab"]),
                               MPP= ifelse(N==0,NA,exp(lnA2)*elastfert*(lab^coef(model)["loglab"])*(N^(elastfert-1))),
                               Npm = (Pn/(constantfactor*Pm))^(1/(elastfert-1)),
                               Ndif = N-Npm)                       

sumzone_olsCD <- db_olsCD %>% group_by(ZONE) %>%
  summarize(
    Ncon=mean(ifelse(N>0, N, NA), na.rm = T),
    N=mean(N, na.rm=T),
    Npm=mean(Npm, na.rm=T),
    MPPmean=mean(MPP[!is.infinite(MPP)], na.rm=T),
    MVC=mean((Pm[!is.infinite(MPP)]*MPP[!is.infinite(MPP)])/Pn[!is.infinite(MPP)], na.rm=T),
    Ndif=mean(Ndif, na.rm=T),
    Number=n())


# Select model
model <- sfaCD

# Note that MPP cannot be calculated for plots with N=0 and are therefore set to 0.
db_sfaCD <- db1 %>% mutate(elastfert  = coef(model)["logN"],
                           imprdf = (coef(model)["impr"]*impr),
                           lnA = coef(model)["(Intercept)"] +
                             (coef(model)["noN"]*noN) +
                             (coef(model)["logarea"]*logarea) +
                             (coef(model)["slope"]*slope) +
                             (coef(model)["elevation"]*elevation) +
                             (coef(model)["rain_wq"]* rain_wq) +
                             (coef(model)["crop_count2"]*crop_count2),
                           lnA2 = lnA + imprdf,
                           constantfactor = exp(lnA2)*elastfert*(lab^coef(model)["loglab"]),
                           MPP= ifelse(N==0,NA,exp(lnA2)*elastfert*(lab^coef(model)["loglab"])*(N^(elastfert-1))),
                           Npm = (Pn/(constantfactor*Pm))^(1/(elastfert-1)),
                           Ndif = N-Npm)                       

sumzone_sfaCD <- db_sfaCD %>% group_by(ZONE) %>%
  summarize(
    Ncon=mean(ifelse(N>0, N, NA), na.rm = T),
    N=mean(N, na.rm=T),
    Npm=mean(Npm, na.rm=T),
    MPPmean=mean(MPP[!is.infinite(MPP)], na.rm=T),
    MVC=mean((Pm[!is.infinite(MPP)]*MPP[!is.infinite(MPP)])/Pn[!is.infinite(MPP)], na.rm=T),
    Ndif=mean(Ndif, na.rm=T),
    Number=n())

# CD case makes sense. Higher MPP in frontier case

# Translog
# core translog inputs
translog_inputs <- c("logN", "loglab")


other_inputs <- c("noN", "slope", "logarea",
                  "elevation", "rain_wq", "impr",
                  "crop_count2")

TL_form_basic <- translog_form("logyld", translog_inputs)
TL_form_full <- paste(TL_form_basic,
                      paste(other_inputs, collapse=" + "), sep=" + ")

# compare ols CD and TL => looks similar
olsTL <- lm(TL_form_full, data=db1)
summary(olsTL)
library(stargazer)
stargazer(olsCD, olsTL, type = "text")

# sfa TL => looks similar
sfaTL <- sfa(TL_form_full, data=db1)
summary(sfaTL)
summary(olsTL)


# MPP and optimal fertilizer use

model <- sfaTL
coef(model)

# MPP sfaTL
db_sfaTL <- db1 %>% mutate(part1 = coef(model)["logN"] +
                             2*coef(model)["I(0.5 * logN^2)"]*logN +
                             coef(model)["logN:loglab"]*loglab,
                          part2 = exp(logyld)/exp(logN),
                          MPP = part1*part2)
                             
sumzone_sfaTL <- db_sfaTL %>% group_by(ZONE) %>%
  summarize(
    Ncon=mean(ifelse(N>0, N, NA), na.rm = T),
    N=mean(N, na.rm=T),
    #Npm=mean(Npm, na.rm=T),
    MPPmean=mean(ifelse(N>0, MPP, NA), na.rm = T),
    Number=n())



# olsTL
model <- olsTL
db_olsTL <- db1 %>% mutate(part1 = coef(model)["logN"] +
                             2*coef(model)["I(0.5 * logN^2)"]*logN +
                             coef(model)["logN:loglab"]*loglab,
                           part2 = exp(logyld)/exp(logN),
                           MPP = part1*part2)

sumzone_olsTL <- db_olsTL %>% group_by(ZONE) %>%
  summarize(
    Ncon=mean(ifelse(N>0, N, NA), na.rm = T),
    N=mean(N, na.rm=T),
    #Npm=mean(Npm, na.rm=T),
    MPPmean=mean(ifelse(N>0, MPP, NA), na.rm = T),
    Number=n())

# Calculate optimal values


# helper function for creating sum

Y_sum <- function(modl){
  
  # get the form of the frontier equation
  form <- strsplit(as.character(formula(modl)), "~ ")[[1]][2]
  
  # combine coefs with variables
  form <- unlist(strsplit(form, " + ", fixed=TRUE))
  coef <- coef(modl)[2:(length(coef(modl))-2)]
  names(coef) <- gsub(":", "*", names(coef))
  coef <- coef[form]
  
  # paste together
  form <- paste(paste(coef,
                      form, sep = " * "), collapse=" + ")
  paste(coef(modl)[1], form, sep=" + ")
}

ysum <- Y_sum(olsTL)
logY <- with(db1, eval(parse(text=ysum)))
Y <- exp(logY)
ysum <- gsub("logN", "log(N)", ysum)

# Function to calculate MPP for TL
# Assumes average values for loglab, Y and relprice, while in fact they vary...
# Possibly better to estimate optimum for every observation and then take average of MPP like before
# Use lapply over plots to do this.


Nopt_f <- function(N, modl){
  ysum <- Y_sum(modl)
  logY <- with(db1, eval(parse(text=ysum)))
  Y <- mean(exp(logY), na.rm=T)
  loglab <- mean(db1$loglab)
  relprice <- mean(db1$relprice)
  MPP = ((coef(modl)["logN"] + 
          2*coef(modl)["I(0.5 * logN^2)"]*log(N) +
          coef(model)["logN:loglab"]*loglab)*
          (Y/N)) -  relprice
  return(MPP)
}


uniroot(function(x){Nopt_f(x, sfaTL)}, interval=c(10, 1000))

N <- seq(1, 1000, length.out=100)
plot(N, Nopt_f(N, sfaTL))
abline(h=0, v = 94)

##############################  

# MPP
MPP_f <- function(N){
  logY <- with(row, eval(parse(text=ysum)))
  Y <- exp(logY)
  logarea <- row$logarea
  loglab <- row$loglab
  logP <- row$logP
  MPP = (coef(modl)[2] + coef(modl)[5] * log(N) + coef(modl)[8] * logarea + coef(modl)[9] * loglab) * (Y/N)
  MPP
}

uniroot(MPP_f, interval=c(10, 1000))


