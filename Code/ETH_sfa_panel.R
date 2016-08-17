#######################################
###### ANALYSIS of ETH panel data #####
#######################################

# CHECK:
# Imputation of dummy variables.
# Mean scale variables (see Henningsen)
# SFA from other package
# Panel estimator
# Including CRE
# relation between residuals and error computed below
# Translog estimation
# Do grubbs test for outliers on N?
# Add texture to soil variables

# From the BID: Addis Ababa, Amhara, Oromiya, SNNP, Tigray, and "other regions" (including Dire Dawa). 
# In those regions with the greatest urban populations, Addis Ababa and Oromiya, 20 EAs were selected; while in all other strata, 15 EAs were selected.

#######################################
############## PACKAGES ETC ###########
#######################################

library(plyr)
library(dplyr)
library(stargazer)
library(broom)
library(DescTools)
library(ggplot2)
library(xtable)
library(frontier)
library(moments)
library(tidyr)
library(openxlsx)
library(frontier)
library(moments)
library(AER)

wdPath <- "D:\\Data\\Projects\\ETHYG"
setwd(wdPath)

source("Code/winsor.r")
options(scipen=999)

#######################################
############## LOAD DATA ##############
#######################################

# Load pooled data
dbP <- readRDS("Cache/Pooled_ETH.rds")

#######################################
############## CLEANING ###############
#######################################


# Select maize plots and head of household
dbP <- filter(dbP, status %in% "HEAD", crop_code %in% 2)

# # Create rel_harv_area variable
# dbP$area_farmer[dbP$area_farmer %in% 0] <- NA
# dbP$harv_area[dbP$harv_area %in% 0] <- NA
# 
# dbP <- dbP %>%
#   mutate( sh_harv_area = harv_area/area_farmer,
#           sh_harv_area = ifelse(sh_harv_area >1, 1, sh_harv_area), # for some farmers the harvested area > plot size. Set to 1
#           rel_harv_area = sh_harv_area * area_gps)

# Create id for plots
dbP <- dbP %>% 
  mutate(id=1:dim(.)[1]) 

# Cleaning and analysis depends strongly on which measure is chosen for area, which is the denominator for many variables.
# there are three possible yield variables. 
# that can be created for the last two waves of data. 
# 1. yld1: above uses the full gps areas as denominator
# 2. yld2: uses harvested area as denominator
# 3. yld3: Uses relative harvest area to correct gps area
# To simplify the code we set these values in this part. Subsequent analysis code can then be used for any definition of yield.

# ETH does not present information  on area harvested for 2011. We use area_gps here yld1)

dbP <- dbP %>% 
  mutate(
    area = area_gps, 
    #area = area_gps,
    yld = crop_qty_harv/area,
    N = N/area)


# As we focus on small scale farmers we restrict area size
dbP <- filter(dbP, area_gps <=10)

# cap yield at 18593 kg/ha, the highest potential yield in ETH (not water limited)
dbP <- filter(dbP, yld <= 18592.85714)

# Restrict attention to plots that use N < 700. 700kg/hectare  represents an upper bound limit associated with inorganic fertilizer use in the United States under irrigated corn conditions (Sheahan & Barett 2014) 
dbP <- filter(dbP, N < 700)

# Filter out plots with zero labour
# NOTE we only have harvest labour
# CHECK might replace this with adult equivalent variable!
dbP <- dbP %>%
  mutate(lab = harv_lab + harv_lab_hire) %>%
  filter(lab >0) %>%
  select(-lab)

# Select relevant variables and complete cases
db0 <- dbP %>% 
  dplyr::select(hhid, ea_id, id, ZONE = REGNAME, REGNAME = ZONENAME, WOREDACODE, KEBELECODE, parcel_id, field_id, holder_id, # ZONE AND REGNAMES reversed to remain consistent with other LSMS
                AEZ, fs,
                SOC, SOC2, ph, ph2, RootDepth, 
                rain_year, rain_wq, 
                #SPEI,
                #YA, YW, YP,
                slope, elevation,
                #nutr_av,
                yld, 
                crop_qty_harv, 
                #sold_qty_kg, sold_qty_gr,
                harv_lab, harv_lab_hire ,
                oxen,
                #ae,
                impr, 
                #fung, herb, # Many missing, not useful
                N, 
                #P, 
                #manure, compost, other_org, # Many missing, not useful
                #crop_stand, cropping,
                #legume, 
                irrig, 
                area, area_tot, area_gps,
                sex, age,
                ed_any, family_size, credit,
                literate, cage, death, N1555,
                dist_hh, dist_road, dist_market, dist_popcenter, #dist_regcap,
                #trans_cost,
                title,
                popEA,
                extension, extension2,
                fert_source,
                #road, cost2small_town, bank, micro_finance, ext_agent,
                crop_count, surveyyear,
                rural, 
                lat, lon)


summary(db0)

#######################################
###### COMPLETE CASES DATABASE ########
#######################################

db0 <- db0 %>%
 do(filter(., complete.cases(.)))


######################################
######## Modify and add variables ####
######################################


# Following Burke
db0$phdum[db0$ph < 55] <- 1
db0$phdum[db0$ph >= 55 & db0$ph <=70] <- 2 # Neutral and best suited for crops
db0$phdum[db0$ph > 70] <- 3
db0$phdum <- factor(db0$phdum)


db0$phdum2[db0$ph2 < 55] <- 1
db0$phdum2[db0$ph2 >= 55 & db0$ph2 <=70] <- 2
db0$phdum2[db0$ph2 > 70] <- 3
db0$phdum2 <- factor(db0$phdum2)

# Recode AEZ into 4 zones
# db0$AEZ2 <- db0$AEZ
# db0$AEZ2 <- mapvalues(db0$AEZ2, from = c("Tropic-warm / semi-arid"), to = c("Tropic-warm"))
# db0$AEZ2 <- mapvalues(db0$AEZ2, from = c("Tropic-warm / sub-humid"), to = c("Tropic-warm"))
# db0$AEZ2 <- factor(db0$AEZ2)

# Crop count > 1
db0$crop_count2[db0$crop_count==1] <- 1
db0$crop_count2[db0$crop_count>1] <- 0

# additional variables
db0 <- db0 %>% mutate (logyld=log(yld),
                       yesN = ifelse(N>0, 1,0), # Dummy when plot does not use fertilizer, following approach of Battese (1997)
                       noN = ifelse(N<=0, 1,0), # Dummy when plot does use fertilizer, following approach of Battese (1997)
                       logN = log(pmax(N, noN)), # maximum of dummy and N following Battese (1997)
                       lab = harv_lab + harv_lab_hire,
                       hirelab_sh = harv_lab_hire/(harv_lab_hire + harv_lab)*100,
                       dumoxen = ifelse(oxen>0, 1,oxen),
                       lab=lab/area,
                       #logae = log(ae),
                       #asset = implmt_value + lvstk2_valu,
                       #assetph=asset/area_tot,
                       #logasset = log(assetph+1),
                       loglab = log(lab+1),
                       logarea = log(area_gps), # area_gps not area because we want to add plot size as proxy for economies of scale
                       rain_wq2 = rain_wq*rain_wq,
                       #pestherb = ifelse(herb==1 | pest==1, 1, 0),
                       #ext = ifelse(ext_dummy_pp==1 | ext_dummy_ph ==1, 1, 0),
                       lograin = log(rain_year),
                       dumfertsource = recode(fert_source, c("'Government' =  1; else = 0")),
                       surveyyear2 = replace(surveyyear==2011, 1, 0))

# Add Translog variables
db0 <- db0 %>% 
  mutate(logN2 = 0.5*logN*logN,
         loglab2 = 0.5*loglab*loglab,
         logNlab = logN*loglab,
         logNimpr = logN*impr,
         loglabimpr = loglab*impr,
         logNirrig = logN*irrig,
         loglabirrig = loglab*irrig,
         logNrain = logN*rain_wq,
         logNoxen = logN*oxen,
         loglaboxen = loglab*oxen)

# Add learning variable: number of fertilizer users in kebele
db0 <- db0 %>%
  group_by(KEBELECODE) %>%
  mutate(fertusers = sum(N >0))
  
# Add CRE variables
# All plot and hh level variables that change over time need to be added.
db0 <- db0 %>%
  group_by(hhid) %>%
  mutate(loglab_bar=mean(loglab, na.rm=TRUE),
         logN_bar=mean(logN, na.rm=TRUE),
         noN_bar=mean(noN, na.rm=TRUE),
         area_bar=mean(area, na.rm=TRUE),
         logarea_bar=mean(logarea, na.rm=TRUE),
         rain_wq_bar=mean(rain_wq, na.rm=TRUE),
         irrig_bar=mean(irrig, na.rm = TRUE),
         impr_bar=mean(impr, na.rm = TRUE),
         oxen_bar=mean(dumoxen, na.rm = TRUE),
         credit_bar=mean(credit,na.rm = TRUE),
         extension_bar=mean(extension, na.rm=TRUE),
         crop_count_bar=mean(crop_count2, na.rm=TRUE),
         slope_bar=mean(slope, na.rm=TRUE),
         family_size_bar=mean(family_size, na.rm=TRUE),
         elevation_bar=mean(slope, na.rm=TRUE)
         ) %>%
  ungroup

  

db0 <- droplevels(db0)
summary(db0)


######################################
##### Get plot specific pricess ######
######################################

# Load and merge price data 
Prices <- readRDS("cache/Prices_ETH.rds")

# Merge with panel data
db1 <- left_join(db0, Prices) %>%
       mutate(relprice = fertilizer/maize) %>%
      rename(Pm = maize, Pn = fertilizer)

# Drop unused levels (e.g. Zanzibar in zone), which are giving problems with sfa
db1 <- droplevels(db1)
saveRDS(db1, "Cache/db1.rds")

############################
###### INPUT DEMAND ########
############################
# CHECK: Add wealth indicator such as Larson, i.e. type of roof or number of livestock.
# Add info on type of farm, mixed, livestock, etc. 
# CHECK OF FARM INCOME

# second stage model - tobit model for adoption
# of nitrogen using HH chars and plot chars as covariates, along
# with CRE time averages to control for unobserved
# household fixed effects
# include soil variables. And other variables leading to higher demand for fertilizer.

N_dem <- tobit(N ~  dist_hh + dist_market +
                    irrig + impr +
                    slope + elevation + 
                    SOC2 + phdum2 +
                    extension + credit +
                    rain_wq + 
                    dumfertsource +
                    sex + age + family_size + literate + ed_any +
                    relprice + 
                    #fertusers +                 
                    surveyyear + crop_count2
                    , data = db1)

N_dem_CRE <- tobit(N ~ dist_hh + dist_market + 
                     irrig +
                     impr +
                     slope + elevation +
                     SOC2 + phdum2 +
                     extension + credit +
                     rain_wq +
                     dumfertsource +
                     loglab + dumoxen + logarea +
                     sex + age + family_size + literate + ed_any + 
                     #log(maize) + log(fertilizer) +
                     relprice + 
                     #fertusers +
                     surveyyear + crop_count2 + 
                     irrig_bar + impr_bar + slope_bar +
                     extension_bar + credit_bar + rain_wq_bar +
                     crop_count_bar + credit_bar + loglab_bar + oxen_bar 
                     , data = db1)



# library(mhurdle)
# NX <- mhurdle(N ~ dist_hh + dist_market + 
#                 irrig +
#                 impr +
#                 slope + elevation +
#                 SOC2 + phdum2 +
#                 extension + credit +
#                 rain_wq +
#                 #rain_wq2+
#                 dumfertsource +
#                 loglab + dumoxen + logarea +
#                 sex + age + family_size + literate + ed_any + 
#                 log(maize) + log(fertilizer) +                  
#                 surveyyear + crop_count2 + 
#                 irrig_bar + impr_bar + crop_count_bar + credit_bar + loglab_bar + oxen_bar |
#                 dist_hh + dist_market + 
#                 irrig +
#                 impr +
#                 slope + elevation +
#                 SOC2 + phdum2 +
#                 extension + credit +
#                 rain_wq +
#                 #rain_wq2+
#                 dumfertsource +
#                 loglab + dumoxen + logarea +
#                 sex + age + family_size + literate + ed_any + 
#                 log(maize) + log(fertilizer) +                  
#                 surveyyear + crop_count2 + 
#                 irrig_bar + impr_bar + crop_count_bar + credit_bar + loglab_bar + oxen_bar 
#               ,dist="ln", data = db1)
# 
# summary(NX)
stargazer(N_dem, N_dem_CRE, type = "text")

summary(N_dem_CRE)

#Take residuals
db1 <- db1 %>%
       mutate(r = residuals(N_dem_CRE, type="deviance"))



#######################################
###### PRODUCTION FUNCTION ############
#######################################

# Cobb Douglas
olsCD <- lm(logyld ~ noN + logN + loglab + 
               dumoxen +
               logarea +
               irrig +
               impr +
               slope + elevation +
               SOC2 + phdum2 + 
               rain_wq + 
               AEZ +
               crop_count2 + surveyyear2 +
               r,
             data = db1)


olsCD_CRE <- lm(logyld ~ noN + logN + loglab + dumoxen + logarea +
               irrig + 
               impr +
               slope + elevation +
               SOC2 + phdum2 + 
               rain_wq + 
               AEZ +
               crop_count2 + surveyyear2 + 
               noN_bar + logN_bar + loglab_bar + logarea_bar + oxen_bar +
               irrig_bar + 
               impr_bar +
               slope_bar + 
               rain_wq_bar +
               crop_count_bar +
               r,
             data = db1)

# Translog function
olsTL2 <- lm(logyld ~ noN + 
               logN + loglab +
               logN2 + loglab2 + logNlab + logNimpr + loglabimpr + logNoxen + loglaboxen +
               #logNirrig +  logNrain + loglabirrig +
               impr +
               dumoxen +
               logarea +
               irrig + 
               slope + elevation +
               SOC2 + phdum2 + 
               rain_wq + 
               AEZ +
               crop_count2 + surveyyear2 + 
               noN_bar + logN_bar + loglab_bar + logarea_bar + oxen_bar +
               irrig_bar + 
               impr_bar +
               slope_bar + elevation_bar +
               crop_count_bar+
               r,
             data = db1)

stargazer(olsCD1, olsCD2, olsTL2, type="text")

# Assess skewness of OLS - should be left skewed which is confirmed.
hist( residuals(olsCD1), 15)
hist( residuals(olsCD2), 15)
skewness(residuals(olsCD1))
skewness(residuals(olsCD2))

# CHECK: NEED TO DEMEAN THE FUNCTION
# Frontier estimation
sfaCD <- sfa(logyld ~ noN + logN +
                loglab +
                dumoxen +
                logarea +
                irrig +
                impr +
                slope + elevation +
                SOC2 + phdum2 +
                rain_wq + 
                AEZ +
                crop_count2 + surveyyear2 +
               r,
              data = db1, maxit = 1500, restartMax = 20, printIter = 1, tol = 0.000001)

summary(sfaCD, extraPar = TRUE)
lrtest(sfaCD)

sfaCD_CRE <- sfa(logyld ~ noN + logN + loglab + logarea + dumoxen +
                irrig + 
                impr +
                slope + elevation +
                SOC2 + phdum2 + 
                rain_wq + 
                AEZ +
                crop_count2 + surveyyear2 + 
                noN_bar + logN_bar + loglab_bar + logarea_bar + oxen_bar +
                irrig_bar + 
                impr_bar +
                slope_bar + 
                #elevation_bar non included because constant over time and plot
                crop_count_bar +
                r,
                data = db1, maxit = 1500, restartMax = 20, tol = 0.000001)

# Translog
sfaTL_CRE <- sfa(logyld ~ noN + logN + loglab + dumoxen +
                   logN2 + loglab2 + logNlab + logNimpr + loglabimpr + logNoxen + loglaboxen + 
                   logarea +
                   irrig + 
                   impr +
                   slope + 
                   SOC2 + phdum2 + 
                   rain_wq + 
                   AEZ +
                   crop_count2 + surveyyear2 + 
                   noN_bar + logN_bar + loglab_bar + logarea_bar + oxen_bar +
                   irrig_bar + 
                   impr_bar +
                   slope_bar + 
                   crop_count_bar +
                   r
                 ,data = db1, maxit = 1500, restartMax = 20, tol = 0.000001)

summary(sfaTL_CRE, extraPar = TRUE)
lrtest(sfaCD_CRE)


# Appears that including FS leads to distance variables to reduce significancy!
sfaCD_CRE_Z <- sfa(logyld ~ noN + logN + loglab + dumoxen + 
                logarea +
                irrig + 
                impr +
                slope + elevation +
                SOC2 + phdum2 + 
                rain_wq + 
                #fs +
                AEZ +
                crop_count2 + surveyyear2 + 
                noN_bar + logN_bar + loglab_bar + logarea_bar + oxen_bar +
                irrig_bar + 
                impr_bar +
                slope_bar +
                rain_wq_bar +
                crop_count_bar + 
                r
              | 
                sex +
                age + 
                title +
                literate +
                ed_any +
                extension +
                credit +
                dist_hh +
                dist_market +
                popEA 
                #hirelab_sh has 38 missing values
                -1
              ,data = db1, maxit = 1500, restartMax = 20, tol = 0.000001)

summary(sfaCD_CRE_Z, extraPar = TRUE)

summary(db1)
lrtest(sfaCD_CRE_Z)


# Compute profit maximizing Pn per zone and other summary statistics

# olsCD_CRE
model <- olsCD_CRE
# Note that MPP cannot be calculated for plots with N=0 and are therefore set to 0.
db_olsCD_CRE <- db1 %>% mutate(elastfert  = coef(model)["logN"],
                               AEZf = ifelse(AEZ=="Tropic - warm / subhumid", coef(model)["AEZTropic - warm / subhumid"], 
                                             ifelse(AEZ=="Tropic - cool / semiarid", coef(model)["AEZTropic - cool / semiarid"],
                                                    ifelse(AEZ=="Tropic - cool / subhumid", coef(model)["AEZTropic - cool / subhumid"],
                                                           ifelse(AEZ=="Tropic - cool / humid", coef(model)["AEZTropic - cool / humid"],
                                                                  0)))),
                               phconstant = ifelse(phdum==2, coef(model)["phdum22"], 
                                                   ifelse(phdum==3, coef(model)["phdum23"],0)),
                               imprdf = (coef(model)["impr"]*impr),
                               lnA = coef(model)["(Intercept)"] +
                                 (coef(model)["noN"]*noN) +
                                 (coef(model)["dumoxen"]*dumoxen) +
                                 (coef(model)["logarea"]*logarea) +
                                 (coef(model)["irrig"]*irrig) +
                                 (coef(model)["slope"]*slope) +
                                 (coef(model)["SOC2"]*SOC2) +
                                 (phconstant) +
                                 (coef(model)["rain_wq"]* rain_wq) +
                                 (coef(model)["crop_count2"]*crop_count2) +
                                 (coef(model)["surveyyear2"]*surveyyear2) +
                                 (coef(model)["noN_bar"]*noN_bar) +
                                 (coef(model)["logN_bar"]*logN_bar) +
                                 (coef(model)["loglab_bar"]*loglab_bar) +
                                 (coef(model)["logarea_bar"]*logarea_bar) +
                                 (coef(model)["oxen_bar"]*oxen_bar) +
                                 (coef(model)["irrig_bar"]*irrig_bar) + 
                                 (coef(model)["impr_bar"]*impr_bar) +
                                 (coef(model)["slope_bar"]*slope_bar) +
                                 (coef(model)["rain_wq_bar"]*rain_wq_bar) +
                                 (coef(model)["crop_count_bar"]*crop_count_bar) +
                                 (coef(model)["r"]*r),
                               lnA2 = lnA + imprdf,
                               constantfactor = exp(lnA2)*elastfert*(lab^coef(model)["loglab"]),
                               MPP= ifelse(N==0,NA,exp(lnA2)*elastfert*(lab^coef(model)["loglab"])*(N^(elastfert-1))),
                               Npm = (Pn/(constantfactor*Pm))^(1/(elastfert-1)),
                               Ndif = N-Npm)                       

sumzone_olsCD_CRE <- db_olsCD_CRE %>% group_by(ZONE) %>%
  summarize(
    Ncon=mean(ifelse(N>0, N, NA), na.rm = T),
    N=mean(N, na.rm=T),
    Npm=mean(Npm, na.rm=T),
    MPPmean=mean(MPP[!is.infinite(MPP)], na.rm=T),
    MVC=mean((Pm[!is.infinite(MPP)]*MPP[!is.infinite(MPP)])/Pn[!is.infinite(MPP)], na.rm=T),
    Ndif=mean(Ndif, na.rm=T),
    Number=n())


# Select model
model <- sfaCD_CRE_Z

# Note that MPP cannot be calculated for plots with N=0 and are therefore set to 0.
db_sfaCD_CRE_Z <- db1 %>% mutate(elastfert  = coef(model)["logN"],
                               AEZf = ifelse(AEZ=="Tropic - warm / subhumid", coef(model)["AEZTropic - warm / subhumid"], 
                                             ifelse(AEZ=="Tropic - cool / semiarid", coef(model)["AEZTropic - cool / semiarid"],
                                                    ifelse(AEZ=="Tropic - cool / subhumid", coef(model)["AEZTropic - cool / subhumid"],
                                                           ifelse(AEZ=="Tropic - cool / humid", coef(model)["AEZTropic - cool / humid"],
                                                                  0)))),
                               phconstant = ifelse(phdum==2, coef(model)["phdum22"], 
                                                   ifelse(phdum==3, coef(model)["phdum23"],0)),
                               imprdf = (coef(model)["impr"]*impr),
                               lnA = coef(model)["(Intercept)"] +
                                 (coef(model)["noN"]*noN) +
                                 (coef(model)["dumoxen"]*dumoxen) +
                                 (coef(model)["logarea"]*logarea) +
                                 (coef(model)["irrig"]*irrig) +
                                 (coef(model)["slope"]*slope) +
                                 (coef(model)["SOC2"]*SOC2) +
                                 (phconstant) +
                                 (coef(model)["rain_wq"]* rain_wq) +
                                 (coef(model)["crop_count2"]*crop_count2) +
                                 (coef(model)["surveyyear2"]*surveyyear2) +
                                 (coef(model)["noN_bar"]*noN_bar) +
                                 (coef(model)["logN_bar"]*logN_bar) +
                                 (coef(model)["loglab_bar"]*loglab_bar) +
                                 (coef(model)["logarea_bar"]*logarea_bar) +
                                 (coef(model)["oxen_bar"]*oxen_bar) +
                                 (coef(model)["irrig_bar"]*irrig_bar) + 
                                 (coef(model)["impr_bar"]*impr_bar) +
                                 (coef(model)["slope_bar"]*slope_bar) +
                                 (coef(model)["rain_wq_bar"]*rain_wq_bar) +
                                 (coef(model)["crop_count_bar"]*crop_count_bar) +
                                 (coef(model)["r"]*r),
                               lnA2 = lnA + imprdf,
                               constantfactor = exp(lnA2)*elastfert*(lab^coef(model)["loglab"]),
                               MPP= ifelse(N==0,NA,exp(lnA2)*elastfert*(lab^coef(model)["loglab"])*(N^(elastfert-1))),
                               Npm = (Pn/(constantfactor*Pm))^(1/(elastfert-1)),
                               Ndif = N-Npm) 

sumzone_sfaCD_CRE_Z<- db_sfaCD_CRE_Z %>% group_by(ZONE) %>%
  summarize(
    Ncon=mean(ifelse(N>0, N, NA), na.rm = T),
    N=mean(N, na.rm=T),
    Npm=mean(Npm, na.rm=T),
    MPPmean=mean(MPP[!is.infinite(MPP)], na.rm=T),
    MVC=mean((Pm[!is.infinite(MPP)]*MPP[!is.infinite(MPP)])/Pn[!is.infinite(MPP)], na.rm=T),
    Ndif=mean(Ndif, na.rm=T),
    Number=n())



##############################
### Calculate yield levels ###
##############################

# set model
model <- sfaCD_CRE_Z

# 1. TEYG: Technical efficient yield gap

# Calculate yield level at 100% TE on the frontier with given inputs
# We estimate the error using the sfa formula and compute Ycor.
# In sfa, Y = TEY + error(v) - efficiency(u) Kumbakar et al. (2015), A practitioner's guide, p.48-49
# We want to filter out/correct for the error. Ycor = TEY - u = Y + e
# Since we do not know the error we calculate Ycor as TEY - u.
# Efficiency as produced by package frontier is defined as TE = exp(-u) so u is -log(TE)
# As the model is in logs, TEY = exp(ln(TEY))
# Ycor = exp(ln(TEY)) - [-log(TE)]

db3 <- db_sfaCD_CRE_Z %>%
          dplyr::select(id, hhid, holder_id, parcel_id, field_id, ZONE, REGNAME, area, crop_count2, lat, lon, lnA, lnA2, noN, yesN, loglab, lab, elastfert, Npm, N, Y=yld) %>%
          mutate(
            Ycor = exp(as.numeric(fitted(model))+log(as.numeric(efficiencies(model)))), 
            err = Ycor-Y,
            TEY = exp(as.numeric(fitted(model))),
            TE = as.numeric(efficiencies(model)),
            resid = as.numeric(resid(model))
            )

# A number of plots have yld higher than the estimated frontier (Y-TEY>0) caused by the random error. 
# A large number of these are plots that do not use fertilizer. They probably have high yield because of measurement error, better soil properties or unknown factors.
above_frontier_check <- filter(db3, Y-TEY>0)
mean(db3$Y-db3$TEY)


# 2. EY: Economic yield 
# Calculate optimal Npm when using Pn
# Note that noN is still part of lnA2 because we assume that plots without N are structurally different from those with N, for instance better soil.
# It is possible that Npm is larger than Npy, which is not possible from a biophysical perspective.
# We cap Npm at Npy.

# Based on experimental plot data (see Excel), we set Npy to 120. CHECK CHANGE FOR ETHIOPIA
Npy <- 120
# Cap Npm
db3 <- mutate(db3, Npm = ifelse(Npm>Npy, Npy, Npm))

db4 <- db3 %>% 
  mutate(EY = exp(
    lnA2 + 
      coef(model)["loglab"]*loglab +
      # coef(model)["logasset"]*logasset + 
      elastfert*log(Npm)))


# 3 PFY: Feasible yield
# To improve this part, we could also argue that: (1) hybrid seeds are used, (2) pestices are used, (3) higher levels of capital and labour are used.
# We assume that all farmers use improved seeds and increase labour with 10%

db5 <- db4 %>%
  mutate(PFY = exp(
    lnA + 
      coef(model)["impr"] +                   # assume all plots use improved seeds
      coef(model)["loglab"]*loglab*1.1 +
      elastfert*log(Npy)))


# 4. PY: Potential yield
# Merge Yield potential with maize plot database
db6 <- dbP %>% dplyr::select(hhid, holder_id, parcel_id, field_id, surveyyear, id, PY = YW) %>% na.omit %>% 
  mutate(PY=PY*1000) %>% left_join(db5,.)

# A large number of plots have missing YW values because region is not covered by GYGA.
# We assume for the moment that country maximum water limited yield (Yw) is reasonable proxy for missing information.
# Might scale this down to see the effect.

GYGA_YW <- 18071.7857142857
db6 <- mutate(db6, PY = ifelse(is.na(PY), GYGA_YW, PY))


#####################################################
### Yield levels consistency check and correction ###
#####################################################

# Because of imputation of TY or measurement error, Yield (Y and Ycor), Technical efficiency yield (TEY), Economic yield (EY) and Feasible yield (UY) 
# can be higher than Potential yield (PYcor). We check for this.

Y_Ycor_check <- filter(db6, Ycor-Y<0)
PY_Y_check <- filter(db6, PY-Y<0)
PY_Y_cor_check <- filter(db6, PY-Ycor<0)
PY_TE_check <- filter(db6, PY-TEY<0)
PY_EY_check <- filter(db6, PY-EY<0)
PY_PFY_check <- filter(db6, PY-PFY<0)

# Compare different yield levels
# Picture shows that PFY is much to high for plots with the lowest PY. This is probably due to the uniform use of Npf of 120 N/ha.
# It would be better to have zone specific Npf values.
ggplot(data = db6, aes(y = PY, x = PFY)) +
  geom_point() +
  #geom_jitter(position=position_jitter(width=.1, height=0))+
  geom_abline(aes(Y = Ycor), slope=1, intercept=0) +
  coord_fixed() +
  scale_y_continuous(limits=c(0, 10000)) +
  scale_x_continuous(limits=c(0, 10000))

# Compare error and resid
# Not clear what resid is? As the following plot shows, TEYG_s = TE :CHECK]
ggplot(data = db6, aes(y = err, x = resid)) +
  geom_point() 

# Compare Sfa TA scores with mannually computed TEYG_s => identical as they should be
db6a <- mutate(db6, TEYG_s =Ycor/TEY)
ggplot(data = db6a, aes(y = TEYG_s, x = TE)) +
  geom_point()

#  We cap all values at PY because we consider this as an absolute potential and recalculate all gaps.
db7 <- mutate(db6, PFY = ifelse(PY-PFY<0, PY, PFY),
              EY = ifelse(PY-EY<0, PY, EY),
              TEY = ifelse(PY-TEY<0, PY, TEY),
              Ycor = ifelse(PY-Ycor<0, PY, Ycor),
              Y = ifelse(PY-Y<0, PY, Y))

#############################
### Yield gap calculation ###
#############################

# Calculate TYG using UY as reference
db8 <- db7 %>% 
  mutate(
    ERROR_l = Ycor-Y,      # Error gap
    ERROR_s = Y/Ycor,      # Error gap
    TEYG_l = TEY-Ycor,     # Technical efficiency yield gap using Ycor as basis
    TEYG_s = Ycor/TEY,     # Technical efficiency yield gap using Ycor as basis
    EYG_l = EY-TEY,        # Economic yield gap
    EYG_s = TEY/EY,        # Economic yield gap
    EUYG_l = PFY-EY,       # Feasible yield gap
    EUYG_s = EY/PFY,       # Feasible yield gap
    TYG_l = PY-PFY,        # Technology yield gap
    TYG_s = PFY/PY,        # Technology yield gap
    YG_l = PY-Y,           # Yield gap
    YG_s = Y/PY,           # Yield gap
    YG_l_Ycor = PY-Ycor,   # Yield gap with Ycor as reference
    YG_s_Ycor = Ycor/PY)   # Yield gap with Ycor as reference

# Consistency check of yield gaps.
# ERROR
ERROR_check <- filter(db8, ERROR_l<0) # Half of observation has a negative error which is what would be expected
mean(db8$ERROR_l)
mean(db8$ERROR_s)

# TEYG
TEYG_check <- filter(db8, TEYG_l<0) # Should be zero
mean(db8$TEYG_s)

# EYG
# A number of plots will have to decrease N use Npm < N. In several cases also plots that do no use N
# will have lower Y when they start using N. This is because there yield can be located above the frontier (based on fertilizer users) because of the positive effect of noN.
# If we believe that these plots are structurally different and do not use fertilizer because of better soils, they will in fact use too much N and have to decrease.
EYG_check <- filter(db8, EYG_l<0)        
mean(db8$EYG_s)

# EUYG
# A number of plots have negative EUYG_l because Npm is larger than Nyw, the nitrogen that is required to achieve Potential yield (Yw).
# We have corrected this so check should be 0.
EUYG_check <- filter(db8, EUYG_l<0)        
mean(db8$EUYG_s)

# TYG
TYG_check <- filter(db8, TYG_l<0)        
mean(db8$TYG_s)

#YG
YG_check <- filter(db8, YG_l<0)        
YG_check2 <- filter(db8, YG_l_Ycor<0)        

# Check if separate yield gaps add up to total yield gap
Overall_check <- db8 %>%
  mutate(check_l = YG_l/(ERROR_l + TEYG_l + EYG_l + EUYG_l + TYG_l), # Note that for a small number of observatios YG_l=0 resulting in 0/0 which is NaN
         check_s = YG_s/(ERROR_s * TEYG_s * EYG_s * EUYG_s * TYG_s),
         check_l2 = YG_l_Ycor/(TEYG_l + EYG_l + EUYG_l + TYG_l),
         check_s2 = YG_s_Ycor/(TEYG_s * EYG_s * EUYG_s * TYG_s))
summary(Overall_check)


# Create database with relevant variables for further analysis
db9 <- dplyr::select(db8, hhid, holder_id, parcel_id, field_id, ZONE, REGNAME, lat, lon, crop_count2, area, Npm, yesN, Y, N, Ycor, TEY, EY, PFY, PY, ERROR_l, ERROR_s, TEYG_l, TEYG_s, EYG_l, EYG_s, 
                     EUYG_l, EUYG_s, TYG_l, TYG_s, YG_l, YG_s, YG_l_Ycor, YG_s_Ycor)

saveRDS(db9, "Cache/db9.rds")

