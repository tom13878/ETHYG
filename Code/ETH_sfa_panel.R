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

# restrict attention to plots that use N < 700. 700kg/hectare  represents an upper bound limit associated with inorganic fertilizer use in the United States under irrigated corn conditions (Sheahan & Barett 2014) 
dbP <- filter(dbP, N < 700)

# Select relevant variables and complete cases
db0 <- dbP %>% 
  dplyr::select(hhid, ea_id, ZONE = REGNAME, REGNAME = ZONENAME, parcel_id, field_id, holder_id, # ZONE AND REGNAMES reversed to remain consistent with other LSMS
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

# Add CRE variables
db0 <- db0 %>%
  group_by(hhid) %>%
  mutate(loglab_bar=mean(loglab, na.rm=TRUE),
         #logae_bar=mean(loglab, na.rm=TRUE),
         logN_bar=mean(logN, na.rm=TRUE),
         noN_bar=mean(noN, na.rm=TRUE),
         area_bar=mean(area, na.rm=TRUE),
         logarea_bar=mean(logarea, na.rm=TRUE),
         irrig_bar=mean(irrig, na.rm = TRUE),
         #legume_bar=mean(legume, na.rm = TRUE),
         #fung_bar=mean(fung, na.rm = TRUE),
         #herb_bar=mean(herb, na.rm = TRUE),
         impr_bar=mean(impr, na.rm = TRUE),
         oxen_bar=mean(dumoxen, na.rm = TRUE),
         credit_bar=mean(credit,na.rm = TRUE),
         crop_count_bar=mean(crop_count2, na.rm=TRUE))

db0 <- droplevels(db0)
summary(db0)


#######################################
############## ADD PRICE DATA #########
#######################################

# Load and merge price data
Prices <- readRDS("./Analysis/ETH/Data/Prices_ETH.rds")
Prices <- Prices %>% dplyr:: select(region_lsms, district, lat, lon, price, type) %>%
  spread(type, price)

# Merge with panel data
db1 <- left_join(db0, Prices)

# Drop unused levels (e.g. Zanzibar in zone), which are giving problems with sfa
db1 <- droplevels(db1)


db0 <- db0 %>%
  do(filter(., complete.cases(.)))

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
                    #rain_wq2
                    sex + age + family_size + literate + ed_any +
                    surveyyear + crop_count2 + dumfertsource
                    , data = db0)

N_dem_CRE <- tobit(N ~ dist_hh + dist_market + 
                     irrig +
                     impr +
                     slope + elevation +
                     SOC2 + phdum2 +
                     extension + credit +
                     rain_wq +
                     #rain_wq2+
                     dumfertsource +
                     loglab + dumoxen + logarea +
                     sex + age + family_size + literate + ed_any + 
                     surveyyear + crop_count2 +
                     irrig_bar + impr_bar + crop_count_bar + credit_bar + loglab_bar + oxen_bar 
                     AEZ
                     , data = db0)

# library(censReg)
# NX <- censReg(N ~ dist_hh + dist_market + 
#                 irrig +
#                 impr +
#                 slope + elevation +
#                 SOC2 + phdum2 +
#                 extension + credit +
#                 rain_wq +
#                 #rain_wq2+
#                 sex + age + family_size + literate + ed_any + fert_source +
#                 surveyyear + crop_count2 +
#                 irrig_bar + impr_bar + crop_count_bar + credit_bar
#               , data = db0)
stargazer(N_dem, N_dem_CRE, type = "text")

summary(N_dem)
#Take residuals
r <- residuals(N_dem_CRE, type="deviance")

require(censReg)

#######################################
###### PRODUCTION FUNCTION ############
#######################################


# Cobb Douglas
olsCD1 <- lm(logyld ~ noN + logN + loglab + 
               dumoxen +
               logarea +
               irrig +
               impr +
               slope + elevation +
               SOC2 + phdum2 + 
               rain_wq + rain_wq2+
               AEZ +
               crop_count2 + surveyyear2,
             data = db0)


olsCD2 <- lm(logyld ~ noN + logN + loglab +
               dumoxen +
               logarea +
               irrig + 
               impr +
               slope + elevation +
               SOC2 + phdum2 + 
               rain_wq + rain_wq2+
               AEZ +
               crop_count2 + surveyyear2 + 
               noN_bar + logN_bar + loglab_bar + oxen_bar +
               irrig_bar + 
               impr_bar +
               crop_count_bar +
               r,
             data = db0)

stargazer(olsCD1, olsCD2, type="text")

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
                rain_wq + rain_wq2+
                AEZ +
                crop_count2 + surveyyear2,
              data = db0, maxit = 1500, restartMax = 20, printIter = 1, tol = 0.000001)

summary(sfaCD, extraPar = TRUE)
lrtest(sfaCD)

sfaCD_CRE <- sfa(logyld ~ noN + logN + 
                loglab +
                dumoxen +
                logarea +
                irrig + 
                impr +
                slope + elevation +
                SOC2 + phdum2 + 
                rain_wq + rain_wq2+
                AEZ +
                crop_count2 + surveyyear2 + 
                noN_bar + logN_bar + loglab_bar +
                irrig_bar + 
                impr_bar +
                crop_count_bar +
                r,
                data = db0, maxit = 1500, restartMax = 20, tol = 0.000001)

summary(sfaCD_CRE, extraPar = TRUE)
lrtest(sfaCD_CRE)

check <- select(db0, sex, literate, cage, ed_any, N1555, family_size, death, dist_hh, dist_road, dist_market, dist_popcenter, dist_regcap,
title, popEA,
road, cost2small_town, bank, micro_finance, ext_agent, extension)
summary(check)

db0X <- db0 %>% mutate(title = as.numeric(ifelse(title == "YES", 1, ifelse(title == "NO", 0, title))),
                       #literate = as.numeric(ifelse(literate == "YES", 1, ifelse(literate == "NO", 0, literate))),
                       ed_any = as.numeric(ifelse(ed_any == "YES", 1, ifelse(ed_any == "NO", 0, ed_any))),
                       #N1555 = as.factor(N1555),
                       extension = as.numeric(ifelse(extension == "YES", 1, ifelse(extension == "NO", 0, extension)))
                       #micro_finance = as.numeric(ifelse(micro_finance == "YES", 1, ifelse(micro_finance == "NO", 0, micro_finance))),
                       #bank = as.numeric(ifelse(bank == "YES", 1, ifelse(bank == "NO", 0, bank)))
                       #popEA = log(popEA+1)
                       )


table(db0X$extension)

check2 <- select(db0X, sex, literate, cage, ed_any, N1555, family_size, death, dist_hh, dist_road, dist_market, dist_popcenter, dist_regcap,
                title, popEA,
                road, cost2small_town, bank, micro_finance, ext_agent, extension)


# Appears that including FS leads to distance variables to reduce significancy!
sfaCD_CRE_Z <- sfa(logyld ~ noN + logN + 
                loglab +
                dumoxen +
                logarea +
                irrig + 
                impr +
                slope + elevation +
                SOC2 + phdum2 + 
                rain_wq + rain_wq2+
                #fs +
                AEZ +
                crop_count2 + surveyyear2 + 
                noN_bar + logN_bar + loglab_bar +
                irrig_bar + 
                impr_bar +
                crop_count_bar + r
              | 
                sex +
                age + 
                title +
                #literate +
                ed_any +
                extension +
                credit +
                # #N1555
                dist_hh +
                dist_market +
                popEA +
                hirelab_sh 
                #micro_finance +
                #cost2small_town 
                #trans_cost
                 # hirelab_sh +
                
                
                # dist_popcenter 
                # #fs +
                 -1
              ,
              data = db0, maxit = 1500, restartMax = 20, tol = 0.000001)

summary(sfaCD_CRE_Z, extraPar = TRUE)
lrtest(sfaCD3)
summary(db0X$popEA)

# Compute profit maximizing Pn per zone and other summary statistics
# Select model
model <- olsCD1

# Note that MPP cannot be calculated for plots with N=0 and are therefore set to 0.
db2_olsCD1 <- db1 %>% mutate(elastfert = ifelse(AEZ2=="Tropic-warm", coef(model)["AEZ2Tropic-warm:logN"], 
                                             ifelse(AEZ2=="Tropic-cool / semi-arid", coef(model)["AEZ2Tropic-cool / semi-arid:logN"],
                                                  ifelse(AEZ2=="Tropic-cool / sub-humid", coef(model)["AEZ2Tropic-cool / sub-humid:logN"], 
                                                         coef(model)["AEZ2Tropic-cool:logN"]))),
                          phconstant = ifelse(phdum==2, coef(model)["phdum22"], 
                                              ifelse(phdum==3, coef(model)["phdum23"],0)),
                          lnA = coef(model)["(Intercept)"] +
                              (coef(model)["noN"]*noN) +
                              (coef(model)["logarea"]*logarea) +
                              (coef(model)["hybrd"]*hybrd) +
                              (coef(model)["manure"]*manure) +
                              (coef(model)["herb"]*herb) +
                              (coef(model)["fung"]*fung) +
                              (coef(model)["legume"]*legume) +
                              (coef(model)["irrig"]*irrig) +
                              (coef(model)["slope"]*slope) +
                              (coef(model)["elevation"]*elevation) +
                              (coef(model)["SOC2"]*SOC2) +
                              (phconstant) +
                              (coef(model)["rain_wq"]* rain_wq) +
                              (coef(model)["rain_wq_2"]*rain_wq_2) +
                              (coef(model)["crop_count2"]*crop_count2),
                          lnA2 = lnA,
                          constantfactor = exp(lnA2)*elastfert*(lab^coef(model)["loglab"]),
                          MPP= ifelse(N==0,NA,exp(lnA2)*elastfert*(lab^coef(model)["loglab"])*(N^(elastfert-1))),
                          Npm = (Pn/(constantfactor*Pm))^(1/(elastfert-1)),
                          Ndif = N-Npm)                       


sumzone_olsCD1<- db2_olsCD1 %>% group_by(region_lsms) %>%
                      summarize(
                        Ncon=mean(ifelse(N>0, N, NA), na.rm = T),
                        N=mean(N, na.rm=T),
                        Npm=mean(Npm, na.rm=T),
                        MPPmean=mean(MPP[!is.infinite(MPP)], na.rm=T),
                        MVC=mean((Pm[!is.infinite(MPP)]*MPP[!is.infinite(MPP)])/Pn[!is.infinite(MPP)], na.rm=T),
                        Ndif=mean(Ndif, na.rm=T),
                        Number=n())

# Select model
model <- sfaCD1_E

# Note that MPP cannot be calculated for plots with N=0 and are therefore set to 0.
db2_sfaCD1_E <- db1 %>% mutate(elastfert = ifelse(AEZ2=="Tropic-warm", coef(model)["AEZ2Tropic-warm:logN"], 
                                                ifelse(AEZ2=="Tropic-cool / semi-arid", coef(model)["AEZ2Tropic-cool / semi-arid:logN"],
                                                       ifelse(AEZ2=="Tropic-cool / sub-humid", coef(model)["AEZ2Tropic-cool / sub-humid:logN"], 
                                                              coef(model)["AEZ2Tropic-cool:logN"]))),
                             phconstant = ifelse(phdum==2, coef(model)["phdum22"], 
                                                 ifelse(phdum==3, coef(model)["phdum23"],0)),
                             lnA = coef(model)["(Intercept)"] +
                               (coef(model)["noN"]*noN) +
                               (coef(model)["logarea"]*logarea) +
                               (coef(model)["hybrd"]*hybrd) +
                               (coef(model)["manure"]*manure) +
                               (coef(model)["herb"]*herb) +
                               (coef(model)["fung"]*fung) +
                               (coef(model)["legume"]*legume) +
                               (coef(model)["irrig"]*irrig) +
                               (coef(model)["slope"]*slope) +
                               (coef(model)["elevation"]*elevation) +
                               (coef(model)["SOC2"]*SOC2) +
                               (phconstant) +
                               (coef(model)["rain_wq"]* rain_wq) +
                               (coef(model)["rain_wq_2"]*rain_wq_2) +
                               (coef(model)["crop_count2"]*crop_count2),
                             lnA2 = lnA,
                             constantfactor = exp(lnA2)*elastfert*(lab^coef(model)["loglab"]),
                             MPP= ifelse(N==0,NA,exp(lnA2)*elastfert*(lab^coef(model)["loglab"])*(N^(elastfert-1))),
                             Npm = (Pn/(constantfactor*Pm))^(1/(elastfert-1)),
                             Ndif = N-Npm)                       



sumzone_sfaCD1_E<- db2_sfaCD1_E %>% group_by(region_lsms) %>%
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
model <- sfaCD1_E

# 1. TEYG: Technical efficient yield gap

# Calculate yield level at 100% TE on the frontier with given inputs
# We estimate the error using the sfa formula and compute Ycor.
# In sfa, Y = TEY + error(v) - efficiency(u) Kumbakar et al. (2015), A practitioner's guide, p.48-49
# We want to filter out/correct for the error. Ycor = TEY - u = Y + e
# Since we do not know the error we calculate Ycor as TEY - u.
# Efficiency as produced by package frontier is defined as TE = exp(-u) so u is -log(TE)
# As the model is in logs, TEY = exp(ln(TEY))
# Ycor = exp(ln(TEY)) - [-log(TE)]

db3 <- db2_sfaCD1_E %>%
          dplyr::select(ea_id2, household_id2, holder_id, parcel_id, field_id, region_lsms, area_im, crop_count2, lat, lon, lnA, lnA2, noN, yesN, loglab, lab, elastfert, Npm, N, Y=yld) %>%
          mutate(
            Ycor = exp(as.numeric(fitted(model))+log(as.numeric(efficiencies(model)))), 
            err = Ycor-Y,
            TEY = exp(as.numeric(fitted(model))),
            TE = as.numeric(efficiencies(model)),
            resid = as.numeric(resid(model))
            )

# A number of plots have yld higher than the estimated frontier (Y-TEY>0) caused by the random error. 
# By far most of these are plots that do not use fertilizer. They probably have high yield because of measurement error, better soil properties or unknown factors.
above_frontier_check <- filter(db3, Y-TEY>0)
mean(db3$Y-db3$TEY)


# 2. EY: Economic yield 
# Calculate optimal Npm when using Pn
# Note that noN is still part of lnA2 because we assume that plots without N are structurally different from those with N, for instance better soil.

# It is possible that Npm is larger than Npy, which is not possible from a biophysical perspective.
# We cap Npm at Npy.

# Based on experimental plot data (see Excel), we set Npy to 150. NB CHANGE FOR ETHIOPIA
Npy <- 150
# Cap Npm
db3 <- mutate(db3, Npm = ifelse(Npm>Npy, Npy, Npm))

db4 <- db3 %>% 
  mutate(EY = exp(
    lnA2 + 
      coef(model)["loglab"]*loglab +
      # coef(model)["logasset"]*logasset + 
      elastfert*log(Npm)))


# 3 PFY: Potential farm yield
db5 <- db4 %>%
  mutate(PFY = exp(
    lnA2 + 
      coef(model)["loglab"]*loglab +
      # coef(model)["logasset"]*logasset + 
      elastfert*log(Npy)))

# 4. PY: Potential yield
# Merge Yield potential with maize plot database
db6 <- CS2013 %>% dplyr::select(household_id2, holder_id, parcel_id, field_id, PY = YW) %>% 
  mutate(PY=PY*1000) %>% left_join(db5, .)

# A large number of plots have missing YW values because region is not covered by GYGA.
# We assume for the moment that country maximum water limited yield (Yw) is reasonable proxy for missing information.
# Might scale this down to see the effect.

###################### CHECK PROBABLY NEEDS TO BE INCREASED WITH A FACTOR 10

GYGA_YW <- 1807.17857142857
db6 <- mutate(db6, PY = ifelse(is.na(PY), GYGA_YW, PY))


#####################################################
### Yield levels consistency check and correction ###
#####################################################

# Because of imputation of TY or measurement error, Yield (Y and Ycor), Technical efficiency yield (TEY), Economic yield (EY) and Unexploitable yield (UY) 
# can be higher than Potential yield (PYcor). We check for this.

Y_Ycor_check <- filter(db6, Ycor-Y<0)
PY_Y_check <- filter(db6, PY-Y<0)
PY_Y_cor_check <- filter(db6, PY-Ycor<0)
PY_TE_check <- filter(db6, PY-TEY<0)
PY_EY_check <- filter(db6, PY-EY<0)
PY_PFY_check <- filter(db6, PY-PFY<0)

# Compare different yield levels
# Picture shows that PFY is much to high for plots with the lowest PY. This is probably due to the uniform use of Npf of 150 N/ha.
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
    EUYG_l = PFY-EY,       # Economically unexploitable yield gap
    EUYG_s = EY/PFY,       # Economically unexploitable yield gap
    TYG_l = PY-PFY,        # Technology yield gap
    TYG_s = PFY/PY,         # Technology yield gap
    YG_l = PY-Y,           # Yield gap
    YG_s = Y/PY,           # Yield gap
    YG_l_Ycor = PY-Ycor,   # Yield gao with Ycor as reference
    YG_s_Ycor = Ycor/PY)   # Yield gap with Ycor as reference

# Consistency check of yield gaps.
# ERROR
ERROR_check <- filter(db8, ERROR_l<0) # Half of observation has a negetive error which is what would be expected # CHECK
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
# Not the case here: CHECK
# We have corrected this so check should be 0.
EUYG_check <- filter(db8, EUYG_l<0)        
mean(db8$EUYG_s)

# TYG
TYG_check <- filter(db8, TYG_l<0)        
mean(db8$TYG_s)

#YG
YG_check <- filter(db8, YG_l<0)        
YG_check2 <- filter(db8, YG_l_Ycor<0)        

# Check if separate yield gaps add up to total yield gap #Check_l has a zero value - CHECK
Overall_check <- db8 %>%
  mutate(check_l = YG_l/(ERROR_l + TEYG_l + EYG_l + EUYG_l + TYG_l), # Note that for a small number of observatios YG_l=0 resulting in 0/0 which is NaN
         check_s = YG_s/(ERROR_s * TEYG_s * EYG_s * EUYG_s * TYG_s),
         check_l2 = YG_l_Ycor/(TEYG_l + EYG_l + EUYG_l + TYG_l),
         check_s2 = YG_s_Ycor/(TEYG_s * EYG_s * EUYG_s * TYG_s))
summary(Overall_check)

# Create database with relevant variables for further analysis
db9 <- dplyr::select(db8, household_id2, holder_id, parcel_id, field_id, region_lsms, lat, lon, crop_count2, area_im, Npm, yesN, Y, N, Ycor, TEY, EY, PFY, PY, ERROR_l, ERROR_s, TEYG_l, TEYG_s, EYG_l, EYG_s, 
                     EUYG_l, EUYG_s, TYG_l, TYG_s, YG_l, YG_s, YG_l_Ycor, YG_s_Ycor)



#######################################
####### FIGURES AND TABLES ############
#######################################

# Number of households and plots per year
HH <- db9 %>% 
  summarize(HH = length(unique(household_id2)),
            Plots = n())


dbsum <- db2_sfaCD1_E %>% dplyr::select(yld, lab, yesN, N, area_im, hybrd, manure, herb, fung, legume, SOC2, rain_wq, slope, elevation, crop_count2, AEZ2)
sumstat <- stargazer(as.data.frame(dbsum), type = "text", digits=2, out=".\\Analysis\\ETH\\Graphs\\summaryStat.html")


# Table with key information per region_lsms and subtotals
Yieldsum <- bind_rows(
  db9 %>% 
    dplyr::select(Y, N, yesN, region_lsms, area_im) %>%
    group_by(region_lsms) %>%
    summarize(Yield = mean(Y),
              Yield_w = (sum(Y*area_im)/sum(area_im)),
              NitrogenUser = round(mean(yesN)*100, digits=1),
              Number=n()),
  db9 %>%
    dplyr::select(Y, N, yesN, area_im, region_lsms) %>%
    summarize(region_lsms = "Total",  
              Yield = mean(Y),
              Yield_w = (sum(Y*area_im)/sum(area_im)),
              NitrogenUser = round(mean(yesN)*100, digits=1),
              Number=n())
) %>% arrange(region_lsms)

# Table with key information per region_lsms and subtotals for pure maize plots
Yieldsum_pure <- bind_rows(
  db9 %>%
    filter(crop_count2==1) %>%
    dplyr::select(Y, N, yesN, region_lsms, area_im) %>%
    group_by(region_lsms) %>%
    summarize(Yield_p = mean(Y),
              Yield_w_p = (sum(Y*area_im)/sum(area_im))
    ),
  db9 %>%
    filter(crop_count2==1) %>%
    dplyr::select(Y, N, yesN, region_lsms, area_im) %>%
    summarize(region_lsms = "Total", 
              Yield_p = mean(Y),
              Yield_w_p = (sum(Y*area_im)/sum(area_im))
    )
)


Nitrogensum <- bind_rows(
  db9 %>% 
    dplyr::select(Y, N, yesN, region_lsms) %>%
    filter(yesN ==1) %>%
    group_by(region_lsms) %>%
    summarize(Nitrogen = mean(N)),
  db9 %>%
    dplyr::select(Y, N, yesN, region_lsms) %>%
    filter(yesN ==1) %>%
    summarize(region_lsms= "Total", Nitrogen = mean(N))
) 


Pricessum <- bind_rows(
  Prices %>% 
    dplyr::select(region_lsms, Pn, Pm) %>% 
    do(filter(., complete.cases(.))) %>%
    group_by(region_lsms) %>%
    summarize(
      NitrogenPrice = mean(Pn, na.rm=T),
      MaizePrice = round(mean(Pm, na.rm=T), digits=0)),
  Prices %>% 
    dplyr::select(region_lsms, Pn, Pm) %>% 
    do(filter(., complete.cases(.))) %>%
    summarize(region_lsms = "Total",
              NitrogenPrice = mean(Pn, na.rm=T),
              MaizePrice = round(mean(Pm, na.rm=T), digits=0))
) 

Zonalsum <- left_join(Yieldsum, Nitrogensum) %>%
  left_join(., Yieldsum_pure) %>%
  left_join(., Pricessum) %>%
  dplyr::select(region_lsms, Number, Yield_w, Yield_w_p, NitrogenUser, Nitrogen, NitrogenPrice, MaizePrice) %>%
  arrange(region_lsms)
Zonalsum <- xtable(Zonalsum, digits = c(0,0,0,0,0,0,0,0,0))
print(Zonalsum, type="html", file=".\\Analysis\\ETH\\Graphs\\Zonal.html")


# SFA table
# Easier to cut and paste in Excel
#sfaTable <- as.data.frame(summary(sfaCD2, extraPar = F)$mleParam)
#print(xtable(sfaTable), type="html", file=".\\Analysis\\TZA\\Graphs\\Sfa.html", digits=3)
summary(sfaCD1_E, extraPar = TRUE)
lrtest(sfaCD1_E)

# Essentials for powerpoint
sfatable <- as.data.frame(summary(sfaCD1_E)$mleParam)[,c(1,3)]
# Production function
sfatable_pf <- sfatable[c(1:22),]
# Efficiency function
sfatable_ef <- sfatable[c(23:29),]

# Table with yield levels
YieldLevels <- bind_rows(
  db9 %>% 
    dplyr::select(region_lsms = region_lsms, Y, Ycor, TEY, EY, PFY, PY, area_im) %>%
    group_by(region_lsms) %>%
    summarize(Y =(sum((Y)*area_im)/sum(area_im)),
              Ycor = (sum((Ycor)*area_im)/sum(area_im)),
              TEY = (sum((TEY)*area_im)/sum(area_im)),
              EY = (sum((EY)*area_im)/sum(area_im)),
              PFY = (sum((PFY)*area_im)/sum(area_im)),
              PY = (sum((PY)*area_im)/sum(area_im))
    ),
  db9 %>% 
    dplyr::select(region_lsms = region_lsms, Y, Ycor, TEY, EY, PFY, PY, area_im) %>%
    summarize(region_lsms = "Total", 
              Y =(sum((Y)*area_im)/sum(area_im)),
              Ycor = (sum((Ycor)*area_im)/sum(area_im)),
              TEY = (sum((TEY)*area_im)/sum(area_im)),
              EY = (sum((EY)*area_im)/sum(area_im)),
              PFY = (sum((PFY)*area_im)/sum(area_im)),
              PY = (sum((PY)*area_im)/sum(area_im)))) %>%
  dplyr::select(region_lsms, Y, Ycor, TEY, EY, PFY, PY)


YieldLevels <- xtable(YieldLevels, digits = c(0,0,0,0,0,0,0,0))
print(YieldLevels, type="html", file=".\\Analysis\\ETH\\Graphs\\YieldLevels.html")


# Table with relative yield gap information per zone
# Note that by definition, YG_s computed by weighting individual YG_s values is not the same as multiplication of weighted TEYG_s etc.
# We therefore calculate YG_s as the product of the weighted components.
ZonalYieldGap_s <- bind_rows(
  db9 %>% 
    dplyr::select(region_lsms = region_lsms, ERROR_s, TEYG_s, EYG_s, EUYG_s, TYG_s, YG_s_Ycor, YG_s, area_im) %>%
    group_by(region_lsms) %>%
    summarize(ERROR_s =(sum((ERROR_s)*area_im)/sum(area_im)),
              TEYG_s = (sum((TEYG_s)*area_im)/sum(area_im)),
              EYG_s = (sum((EYG_s)*area_im)/sum(area_im)),
              EUYG_s = (sum((EUYG_s)*area_im)/sum(area_im)),
              TYG_s = (sum((TYG_s)*area_im)/sum(area_im)),
              YG_s = (sum((YG_s)*area_im)/sum(area_im)),
              YG_s_Ycor = (sum((YG_s_Ycor)*area_im)/sum(area_im)),
              ERROR = (1-ERROR_s)*100,
              TEYG = (1-TEYG_s)*100,
              EYG = (1-EYG_s)*100,
              EUYG = (1-EUYG_s)*100,
              TYG = (1-TYG_s)*100,
              YG = (1-(ERROR_s*TEYG_s*EYG_s*EUYG_s*TYG_s))*100,
              YG_Ycor = (1-(TEYG_s*EYG_s*EUYG_s*TYG_s))*100
    ),
  db9 %>% 
    dplyr::select(region_lsms = region_lsms, ERROR_s, TEYG_s, EYG_s, EUYG_s, TYG_s, YG_s_Ycor, YG_s, area_im) %>%
    summarize(region_lsms = "Total", 
              ERROR_s =(sum((ERROR_s)*area_im)/sum(area_im)),
              TEYG_s = (sum((TEYG_s)*area_im)/sum(area_im)),
              EYG_s = (sum((EYG_s)*area_im)/sum(area_im)),
              EUYG_s = (sum((EUYG_s)*area_im)/sum(area_im)),
              TYG_s = (sum((TYG_s)*area_im)/sum(area_im)),
              YG_s = (sum((YG_s)*area_im)/sum(area_im)),
              YG_s_Ycor = (sum((YG_s_Ycor)*area_im)/sum(area_im)),
              ERROR = (1-ERROR_s)*100,
              TEYG = (1-TEYG_s)*100,
              EYG = (1-EYG_s)*100,
              EUYG = (1-EUYG_s)*100,
              TYG = (1-TYG_s)*100,
              YG = (1-(ERROR_s*TEYG_s*EYG_s*EUYG_s*TYG_s))*100,
              YG_Ycor = (1-(TEYG_s*EYG_s*EUYG_s*TYG_s))*100)) %>%
  dplyr::select(region_lsms, TEYG, EYG, EUYG, TYG, YG = YG_Ycor)


ZonalYieldGap_s <- xtable(ZonalYieldGap_s, digits = c(0,0,0,0,0,0,0))
print(ZonalYieldGap_s, type="html", file=".\\Analysis\\ETH\\Graphs\\ZonalYG_s.html")

# Table with absolute yield gap information per zone
# Note that by definition, YG_s computed by weighting individual YG_s values is not the same as multiplication of weighted TEYG_s etc.
# We therefore calculate YG_s as the product of the weighted components.
ZonalYieldGap_l <- bind_rows(
  db9 %>% 
    dplyr::select(region_lsms = region_lsms, ERROR_l, TEYG_l, EYG_l, EUYG_l, TYG_l, YG_l_Ycor, YG_l, area_im) %>%
    group_by(region_lsms) %>%
    summarize(ERROR_l =(sum((ERROR_l)*area_im)/sum(area_im)),
              TEYG_l = (sum((TEYG_l)*area_im)/sum(area_im)),
              EYG_l = (sum((EYG_l)*area_im)/sum(area_im)),
              EUYG_l = (sum((EUYG_l)*area_im)/sum(area_im)),
              TYG_l = (sum((TYG_l)*area_im)/sum(area_im)),
              YG_l = (sum((YG_l)*area_im)/sum(area_im)),
              YG_l_Ycor = (sum((YG_l_Ycor)*area_im)/sum(area_im)),
              YG_lcheck = (ERROR_l+TEYG_l+EYG_l+EUYG_l+TYG_l)),
  db9 %>% 
    dplyr::select(region_lsms = region_lsms, ERROR_l, TEYG_l, EYG_l, EUYG_l, TYG_l, YG_l_Ycor, YG_l, area_im) %>%
    summarize(region_lsms = "Total", 
              ERROR_l =(sum((ERROR_l)*area_im)/sum(area_im)),
              TEYG_l = (sum((TEYG_l)*area_im)/sum(area_im)),
              EYG_l = (sum((EYG_l)*area_im)/sum(area_im)),
              EUYG_l = (sum((EUYG_l)*area_im)/sum(area_im)),
              TYG_l = (sum((TYG_l)*area_im)/sum(area_im)),
              YG_l = (sum((YG_l)*area_im)/sum(area_im)),
              YG_l_Ycor = (sum((YG_l_Ycor)*area_im)/sum(area_im)),
              YG_lcheck = (ERROR_l+TEYG_l+EYG_l+EUYG_l+TYG_l))
) %>%
  #dplyr::select(-ERROR_l, -YG_lcheck, - YG_l)

# # Addtional code to compute percentage of YG
#%>%
   mutate(
        TEYG = 100*TEYG_l/YG_l_Ycor,
        EYG = 100*EYG_l/YG_l_Ycor,
        EUYG = 100*EUYG_l/YG_l_Ycor,
        TYG = 100*TYG_l/YG_l_Ycor,
        YG = 100*(TEYG_l + EYG_l + EUYG_l + TYG_l)/YG_l_Ycor) %>%
  dplyr::select(-ERROR_l:-YG_lcheck)


ZonalYieldGap_l <- xtable(ZonalYieldGap_l, digits = c(0,0,0,0,0,0,0))
print(ZonalYieldGap_l, type="html", file=".\\Analysis\\ETH\\Graphs\\ZonalYG_l.html")

# Calculation of potential increase in production when gap is closed on the basis of sample
GapClose1 <- mutate(db9, PROD = Y * area_im,
                    ERROR_close = ERROR_l*area_im,
                    TEYG_close = TEYG_l*area_im,
                    EYG_close = EYG_l*area_im,
                    EUYG_close = EUYG_l*area_im,
                    TYG_close = TYG_l*area_im,
                    POTPROD = PROD + ERROR_close + TEYG_close + EYG_close + TYG_close + EUYG_close)

# Total increase in yield per year
GapClose1a <- GapClose1 %>%
  #group_by(region_lsms) %>%
  summarize(PROD = sum(PROD, na.rm=T)/1000000,
            ERROR_close = sum(ERROR_close, na.rm=T)/1000000,
            TEYG_close = sum(TEYG_close, na.rm=T)/1000000,
            EYG_close = sum(EYG_close, na.rm=T)/1000000,
            EUYG_close = sum(EUYG_close, na.rm=T)/1000000,
            TYG_close = sum(TYG_close, na.rm=T)/1000000,
            POTPROD = sum(POTPROD, na.rm=T)/1000000)


# Calculation of potential increase in production when gap is closed on the basis of SPAM and FAO data.
# As the yield in the LSMS is much lower than the FAO/SPAM yield and also that of the LSMS we apply the different yield shares as found above to the base yield of SPAM.
# We use Ycor as base. This means we assume there was an error (e) and this is corrected for.

SPAMData <- read.csv("./Analysis/ETH/Data/SPAMData_ETH.csv") %>%
  rename(Y_SPAM = yield, PROD = TargetProduction, region_lsms = zone)
  

# Closing of yield gaps per zone
# Note that for some regions, notably those with very low potential yield (Central and Western), closing TEY and EY already results in
# Closing the gap. To avoid negative closing the EYG, EUYG and TYG are capped.
# The reason for overshooting can be caused by a variety of factors, including mismeasurement. Most likely is that Nyp is too high for regions
# With a very low potential. The all over impact is low as the involved regions have very limited maize production.


GapClose2 <- db9 %>% 
  group_by(region_lsms) %>%
  summarize(
    TEYG_s = sum(TEYG_s*area_im)/sum(area_im),
    EYG_s = sum(EYG_s*area_im)/sum(area_im),
    #TYG_s = (sum((TYG_s)*area_im)/sum(area_im)), # TYG_s based on LSMS yield, not used
    #YG_s = (sum((YG_s)*area_im)/sum(area_im)), # YG_s based on LSMS yield, not used
    EUYG_s = sum(EUYG_s*area_im)/sum(area_im),
    PY = mean(PY, na.rm=T)) %>% # Average of potential yield from GYGA
  left_join(SPAMData, .) %>%
  mutate(
    TEY_SPAM = Y_SPAM/TEYG_s, # TEY using SPAM yield as reference
    EY_SPAM = TEY_SPAM/EYG_s, # EY using SPAM yield as reference
    EY_SPAM = ifelse(PY-EY_SPAM<0, PY, EY_SPAM), # Correct EY if impact of TEYG and EYG results in yield larger than PY.
    EYG_s_SPAM =  TEY_SPAM/EY_SPAM, # Recalculate EYG_s
    UY_SPAM = EY_SPAM/EUYG_s, # UY using SPAM yield as reference
    UY_SPAM = ifelse(PY-UY_SPAM<0, PY, UY_SPAM), # Correct UY if impact of TEYG and EYG results in yield larger than PY.
    EUYG_s_SPAM =  EY_SPAM/UY_SPAM, # Recalculate UEYG_s
    TYG_s_SPAM = UY_SPAM/PY, # Recalculate TYG_s 
    check = TEYG_s*EYG_s_SPAM*EUYG_s_SPAM*TYG_s_SPAM, #check if multiplication of different parts is the same as total
    YG_s = Y_SPAM/PY, # YG_s using SPAM yield as reference
    PTEYG = PROD/TEYG_s, # Total production when TEYG is closed
    PEYG = PTEYG/EYG_s_SPAM, # Total production when EYG is closed
    PEUYG = PEYG/EUYG_s_SPAM, # Total production when EUYG is closed
    PTYG = PEUYG/TYG_s_SPAM, # Total production when TYG is closed
    POTPROD = PROD/YG_s, # Total production when YG is closed
    TEYG_close = PTEYG - PROD, # Additional production when closing TEYG
    EYG_close = PEYG - PTEYG, # Additional production when closing EYG
    EUYG_close = PEUYG - PEYG, # Additional production when closing EUYG
    TYG_close = POTPROD - PEUYG) %>%
  mutate(check2 = TEYG_close + EYG_close + EUYG_close + TYG_close+PROD)

GapClose2a <- GapClose2 %>% 
  summarize(PROD = sum(PROD/1000000, na.rm = T), # in million tons
            TEYG_close = sum(TEYG_close/1000000, na.rm = T),
            EYG_close = sum(EYG_close/1000000, na.rm = T),
            TYG_close = sum(TYG_close/1000000, na.rm = T),
            EUYG_close = sum(EUYG_close/1000000, na.rm = T), 
            POTPROD = sum(POTPROD/1000000, na.rm = T)) %>%
  mutate(check2 = TEYG_close + EYG_close + EUYG_close + TYG_close+PROD)


# http://www.r-bloggers.com/waterfall-plots-in-r/
# Create database for waterfall plot
# Add error, which is very small to target production
wf.df <- GapClose2a %>% 
  dplyr::select(PROD,  TEYG_close, EYG_close, TYG_close, EUYG_close, POTPROD)

wf.df <- as.data.frame(t(wf.df)) %>%
  mutate(category =c("Actual production", " Closing \n technical efficiency \n yield gap", " Closing \n economic \n yield gap",
                     " Closing \n technical \n yield gap", " Closing \n economically unexploitable \n yield gap", "Potential \n production"),
         sector = category) %>%
  rename(value = V1)

# Create waterfall plot
cbPalette <- c("#009E73", "#CC79A7", "#0072B2", "#D55E00", "black", "#999999")
waterfall_f(wf.df)

## Determines the spacing between columns in the waterfall chart
offset <- 0.3

waterfall <- waterfall_f(wf.df, offset=offset) +
  scale_fill_manual(guide="none", values=cbPalette)+
  labs(x="", y="Maize production (million tons)") +
  scale_y_continuous(breaks=seq(0, 30, 5), labels = comma) +
  theme_classic() 


print(waterfall)
library(Cairo)
ggsave(plot = waterfall, ".\\Analysis\\ETH\\Graphs\\Waterfall.png", height = 150, width = 200, type = "cairo-png", units="mm")


# Distribution of relative yield gaps
db11 <- db9%>%
  dplyr::select(region_lsms, ERROR_s, TEYG_s, EYG_s, EUYG_s, TYG_s, YG_s_Ycor) %>%
  gather(yieldgap, value, ERROR_s:YG_s_Ycor) %>%
  filter(yieldgap!="ERROR_s") %>% 
  droplevels() %>%
  mutate(value = (1-value)*100,
         yieldgap = factor(yieldgap, levels = c("TEYG_s", "EYG_s", "TYG_s", "EUYG_s", "YG_s_Ycor")))


# Note that average error=0 and therefore not interesting to show.
# Show plot with minimum hinges and no outliers.
boxplot <-ggplot(data=db11, aes(x=yieldgap, y=value)) +
  #geom_boxplot(fill=c("#D55E00", "#009E73", "#0072B2", "#CC79A7", "#999999"), outlier.colour = NA) +
  geom_boxplot(fill=c("#D55E00", "#009E73", "#0072B2", "#CC79A7", "#999999"), outlier.colour = "black") +
  stat_boxplot(geom ='errorbar') +
  guides(fill=FALSE) +
  stat_summary(fun.y=mean, geom="point", shape=5, size=4) +
  theme_classic() +
  #geom_jitter(position=position_jitter(width=.1, height=0)) +
  labs(x="", y="Yield gap (%)")+ 
  scale_x_discrete(breaks=c("TEYG_s", "EYG_s", "TYG_s", "EUYG_s", "YG_s_Ycor"), 
                   labels=c("Techical efficiency \n yield gap", "Economic \n yield gap", 
                            "Technical \n yield gap", " Economically \n unexploitable \n yield gap", "Yield gap")) 
#facet_wrap(~zone) +
#scale_y_continuous(breaks=seq(0, 100, 25)) 

# Rescale axes.
#sts <- boxplot.stats(db11$value)$stats  # Compute lower and upper whisker limits
#boxplot = boxplot + coord_cartesian(ylim = c(-5,max(sts)*1.05))
boxplot
ggsave(plot = boxplot, ".\\Analysis\\ETH\\Graphs\\Distribution.png", height = 150, width = 200, type = "cairo-png", units="mm")

# Distribution of absolute yield gaps
db12 <- db9%>%
  dplyr::select(region_lsms, ERROR_l, TEYG_l, EYG_l, EUYG_l, TYG_l, YG_l_Ycor) %>%
  gather(yieldgap, value, ERROR_l:YG_l_Ycor) %>%
  filter(yieldgap!="ERROR_l") %>% 
  droplevels() %>%
  mutate(yieldgap = factor(yieldgap, levels = c("TEYG_l", "EYG_l", "TYG_l", "EUYG_l", "YG_l_Ycor")))

# Note that average error=0 and therefore not interesting to show.
# Show plot with minimum hinges and no outliers.
boxplot2 <-ggplot(data=db12, aes(x=yieldgap, y=value)) +
  geom_boxplot(fill=c("#D55E00", "#009E73", "#0072B2", "#CC79A7", "#999999"), outlier.colour = NA) +
  #geom_boxplot(fill=c("#D55E00", "#009E73", "#0072B2", "#CC79A7", "#999999"), outlier.colour = "black") +
  stat_boxplot(geom ='errorbar') +
  guides(fill=FALSE) +
  stat_summary(fun.y=mean, geom="point", shape=5, size=4) +
  theme_classic() +
  labs(x="", y="Yield gap (%)")+ 
  scale_x_discrete(breaks=c("TEYG_l", "EYG_l", "TYG_l", "EUYG_l", "YG_l_Ycor"), 
                   labels=c("Techical efficiency \n yield gap", "Economic \n yield gap", 
                            "Technical \n yield gap", " Economically \n unexploitable \n yield gap", "Yield gap"))
#+ scale_y_continuous(breaks=seq(0, 100, 25)) 

# Rescale axes.
boxplot2
ggsave(plot = boxplot2, ".\\Analysis\\ETH\\Graphs\\Distribution2.png", height = 150, width = 200, type = "cairo-png", units="mm")


# Maps with GYGA yield potential and plot information
# transform shapefile in dataframe for ggplot. rownames are used as ID for the countries. Using ID gives strange results. 
# Additional data is linked back again
library(rgdal)
library(sp)
library(raster)
library(foreign)
library(gdata)

# GYGA DATA
GYGApath <- "D:\\Data\\IPOP\\GYGA\\"

dsn=paste(GYGApath, "\\CZ_SubSaharanAfrica\\CZ_AFRGYGACNTRY.shp", sep="")
ogrListLayers(dsn)
ogrInfo(dsn, layer="CZ_AFRGYGACNTRY")
GYGA.Africa<-readOGR(dsn, layer = "CZ_AFRGYGACNTRY")
projection(GYGA.Africa) # check projection
GYGA.Africa <- spTransform(GYGA.Africa, CRS("+proj=longlat +datum=WGS84"))

# Get GYGA
GYGA.country.yield.data <- read.xls(paste(GYGApath, "GygaEthiopia.xlsx", sep="\\"), sheet=3)

# Select data for maize
GYGA.country.yield.data <- subset(GYGA.country.yield.data, CROP=="Rainfed maize")

# Cut out ETH from GYGA map
GYGA.country <- GYGA.Africa[GYGA.Africa$REG_NAME=="Ethiopia",]

# Link yield gap data
# in order to merge data with spatialpolygondataframe the row.names need to be the same.
# For this reason I first merge the additionald data and then create a spatialpolygondataframe
# again ensuring that the rownames are preserved.

GYGA.country.data <- as(GYGA.country, "data.frame")
GYGA.country.data$id <-row.names(GYGA.country.data)
GYGA.country.data <- merge(GYGA.country.data, GYGA.country.yield.data[,c(1:8)], by.x=c("GRIDCODE"), by.y=c("CLIMATEZONE"), all.x=TRUE, sort=FALSE)
row.names(GYGA.country.data)<-GYGA.country.data$id
GYGA.country <- SpatialPolygonsDataFrame(as(GYGA.country, "SpatialPolygons"),
                                     data=GYGA.country.data)

# Maps with GYGA yield potential and plot information
# transform shapefile in dataframe for ggplot. rownames are used as ID for the countries. Using ID gives strange results. 
# Additional data is linked back again
GYGA.country.fort<- fortify(GYGA.country) 
GYGA.country.fort <- merge(GYGA.country.fort, GYGA.country.data, by="id")
GYGA.country.fort$yieldclass <- cut(GYGA.country.fort$YW, breaks=c(6, 8.5, 11, 13.5, 16, 19))
meanYield <- ddply(db9,.(lon, lat), summarize, meanYield = (sum(Y*area_im)/sum(area_im))/1000)
meanYield$meanYield2 <- cut(meanYield$meanYield, breaks=c(0, 1, 2, 3, 11))
library(ggthemes)
GYGA_LSMS <- ggplot()+
  geom_polygon(data=GYGA.country.fort, aes(x=long, y=lat, group=group, fill=yieldclass), colour="black")+
  geom_polygon(data=subset(GYGA.country.fort, is.na(YW)), aes(x=long, y=lat, group=group), fill="white", colour="black") +
  scale_fill_discrete(name="Potential water\nlimited yield (tons)") +
  geom_point(data=meanYield, aes(x=lon, y=lat, size=(meanYield2)), colour="black")+
  scale_size_manual(name="Average yield (tons)", values=c(1, 2, 3, 4)) +
  coord_equal()+
  labs(x="", y="")+
  theme_classic() +
  theme(legend.key=element_blank(),
        line = element_blank(),
        axis.text = element_blank())

GYGA_LSMS
ggsave(plot = GYGA_LSMS, ".\\Analysis\\ETH\\Graphs\\GYGA_LSMS.png", height = 150, width = 200, type = "cairo-png", units="mm")

# Zonal map with community yield levels
# read in map of tanzania as a SpatialPolygonsDataFrame
library(raster)
countryMap <- getData('GADM', country = "ETH", level = 1) 

# Rename zones using LSMS names
countryMap@data$region_lsms <- countryMap@data$NAME_1
countryMap@data$region_lsms[countryMap@data$NAME_1 %in% c("Oromia")] <- "Oromiya"
countryMap@data$region_lsms[countryMap@data$NAME_1 %in% c("Somali")] <- "Somalie"
countryMap@data$region_lsms[countryMap@data$NAME_1 %in% c("Benshangul-Gumaz")] <- "Benishangul Gumuz"
countryMap@data$region_lsms[countryMap@data$NAME_1 %in% c("Southern Nations, Nationalities and Peoples")] <- "SNNP"
countryMap@data$region_lsms[countryMap@data$NAME_1 %in% c("Gambela Peoples")] <- "Gambella"
countryMap@data$region_lsms[countryMap@data$NAME_1 %in% c("Harari People")] <- "Harari"
countryMap@data$region_lsms <- factor(countryMap@data$region_lsms)

# Remove Addis Ababa and Dire Dawa
countryMap <- countryMap[!(countryMap@data$region_lsms %in% c("Addis Ababa", "Dire Dawa")),]
plot(countryMap)

#  fortify spatial data to use with ggplot and join using join functions from dplyr
#    The join is on id, make sure all ids are character vectors
tf <- fortify(countryMap)
countryMap@data <- rename(countryMap@data, id = ID_1)
countryMap@data$id <- as.character(countryMap@data$id)
tf2 <- left_join(tf, countryMap@data)

# Use ggplot to plot map of Tanzania, specifying the labels and choosing nice colours
#    from the RColorBrewer package
library(RColorBrewer)
#display.brewer.all()

p4 <- ggplot()+
  geom_polygon(data=tf2, aes(x=long, y=lat, group=group, fill=region_lsms), colour="black")+
  geom_point(data=meanYield, aes(x=lon, y=lat, size=(meanYield2)), colour="black")+
  scale_fill_brewer(name = "Zones", palette = "Set1") +
  scale_size_manual(name="Average yield (tons)", values=c(1.5, 2.5, 3.5, 4,5)) +
  coord_equal()+
  labs(x="", y="")+
  theme_classic() +
  theme(legend.key=element_blank(),
      line = element_blank(),
      axis.text = element_blank())
p4 

ggsave(plot = p4, ".\\Analysis\\ETH\\Graphs\\datamap.png", height = 150, width = 200, type = "cairo-png", units="mm")

save.image(file=".\\Analysis\\ETH\\Report\\IMAGINE_GHA\\ETH_ws.RData")
  

