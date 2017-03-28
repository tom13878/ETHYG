#'=================================================================================================
#' Project:  IMAGINE ETH
#' Subject:  Clean Ethiopia data for analysis
#' Author:   Michiel van Dijk & Tom Morley
#' Contact:  michiel.vandijk@wur.nl, Tomas.morley@wur.nl
#' Output:   Cleaned dataset for analysis
#'=================================================================================================

library(rprojroot)
library(dplyr)

# set working directory for cached files
root <- find_root(is_rstudio_project)

# Load 2013 data
source(file.path(root, "Code/ETH_2013.r"))

# prices won't join properly without
# this because prices use year 1 hhid
ETH2013$household_id <- ifelse(is.na(ETH2013$household_id), ETH2013$household_id2, ETH2013$household_id)
ETH2013$individual_id <- ifelse(is.na(ETH2013$individual_id), ETH2013$individual_id2, ETH2013$individual_id)
ETH2013$ea_id <- ifelse(is.na(ETH2013$ea_id), ETH2013$ea_id2, ETH2013$ea_id)

# continue with dbP
dbP <- ETH2013
dbP <- select(dbP, hhid=household_id, indidy=individual_id, everything())

# Select maize plots and head of household
dbP <- filter(dbP, status %in% "HEAD", crop_code %in% 2)

# Create id for plots
dbP <- dbP %>% 
  mutate(id=1:dim(.)[1]) 

# create a yield and nitrogen application variable
dbP <- dbP %>% 
  mutate(
    area = area_gps, 
    yld = crop_qty_harv/area,
    N = N/area,
    seedha = seed_q/area)


# As we focus on small scale farmers we restrict area size
dbP <- filter(dbP, area_gps <=10)

# cap yield at 18593 kg/ha, the highest potential yield in ETH (not water limited)
dbP <- filter(dbP, yld <= 18592.85714)

# restrict attention to plots that use N < 700. 700kg/hectare
# represents an upper bound limit associated with inorganic
# fertilizer use in the United States under irrigated corn
# conditions (Sheahan & Barett 2014) 
dbP <- filter(dbP, N < 700)

# Select relevant variables and complete cases
db0 <- dbP %>%
  dplyr::select(hhid, ea_id, ZONE = REGNAME, REGNAME = ZONENAME, parcel_id, field_id, holder_id,
                AEZ, fs,
                rain_year, rain_wq,
                ph, ph2, SOC, SOC2, rain_CRU, RootDepth,
                fs, YA, YW, YP, AI, TS, GGD, twi,
                slope, elevation,
                nutr_av,
                yld,
                harv_lab, harv_lab_hire ,
                impr,
                fung, herb,
                N, P, seedha,
                off_farm_income, 
                manure, compost, other_org,
                crop_stand, cropping,
                legume, irrig,
                area, area_tot, area_gps,
                sex, age,
                literate, cage, ed_any, N1555, family_size, death,
                dist_hh, dist_road, dist_market, dist_popcenter, dist_regcap,
                title, credit, rotation, extension2, popEA, 
                popEA, HHEA, oxen, fallow10, fallow_year,
                road, cost2small_town, cost2large_town, bank, micro_finance, ext_agent, extension,
                crop_count, surveyyear,
                rural,
                lat, lon)

# Dummy for oxen
db0$dumoxen <- ifelse(db0$oxen > 0, 1, 0)

#db0$phduml55 <- ifelse(db0$ph < 55, 1, 0)
db0$phdum55_2_70 <- ifelse(db0$ph >= 55 & db0$ph <=70, 1, 0) # Neutral and best suited for crops
db0$phdum_gt70 <- ifelse(db0$ph > 70, 1, 0)  
#db0$phdum <- factor(db0$phdum)


# db0$phdum2[db0$ph2 < 55] <- 1
# db0$phdum2[db0$ph2 >= 55 & db0$ph2 <=70] <- 2
# db0$phdum2[db0$ph2 > 70] <- 3
# db0$phdum2 <- factor(db0$phdum2)

# Crop count > 1
db0$crop_count2[db0$crop_count==1] <- 1
db0$crop_count2[db0$crop_count>1] <- 0

# additional variables
db0 <- db0 %>% mutate (logyld=log(yld),
                       yesN = ifelse(N>0, 1,0), # Dummy when plot does not use fertilizer, following approach of Battese (1997)
                       noN = ifelse(N<=0, 1,0), # Dummy when plot does use fertilizer, following approach of Battese (1997)
                       noP = ifelse(P<=0, 1,0), # Dummy when plot does use fertilizer, following approach of Battese (1997)
                       logN = log(pmax(N, noN)), # maximum of dummy and N following Battese (1997)
                       logP = log(pmax(P, noP)), # maximum of dummy and N following Battese (1997)
                       lab = harv_lab + harv_lab_hire,
                       logslope = log(slope+1), # log transformation is more suitable for slope
                       hirelab_sh = harv_lab_hire/(harv_lab_hire + harv_lab)*100,
                       lab=lab/area,
                       loglab = log(lab+1),
                       logseed = log(seedha + 1),
                       logarea = log(area_gps), # area_gps not area because we want to add plot size as proxy for economies of scale
                       rain_wq2 = rain_wq*rain_wq,
                       rain_year2 = rain_year*rain_year,
                       lograin = log(rain_year),
                       #sex = as.numeric(ifelse(sex == "MALE", 0, 1)),
                       surveyyear2 = replace(surveyyear==2011, 1, 0))

db0 <- droplevels(db0)

# Load price data 
Prices <- readRDS(file.path(root, "cache/Prices_ETH.rds"))
#%>% select(hhid, fertilizer, maize)

# Merge with panel data
db1 <- left_join(db0, Prices) %>%
  mutate(relprice = fertilizer/maize) %>%
  rename(Pm = maize, Pn = fertilizer)

# Drop unused levels (e.g. Zanzibar in zone), which are giving problems with sfa
db1 <- droplevels(db1)

# the crop stand variable has problems
db1$crop_stand <- ifelse(db1$crop_stand %in% c("NAN", "OTHER LAND USE (SPECIFY)"), NA, db1$crop_stand)
db1$crop_stand <- ifelse(db1$crop_stand %in% c("PURESTAND"), "PURE STAND", db1$crop_stand)
db1$crop_stand <- ifelse(db1$crop_stand %in% c("PURE STAND"), 1,
                         ifelse(db1$crop_stand %in% c("MIXED CROP"), 0, NA))
db1$micro_finance <- ifelse(db1$micro_finance %in% c("YES"), 1,
                         ifelse(db1$micro_finance %in% c("NO"), 0, NA))
db1$manure <- ifelse(db1$manure %in% 2, 1,
                            ifelse(db1$manure %in% 1, 0, NA))
db1$compost <- ifelse(db1$compost %in% 2, 1,
                     ifelse(db1$compost %in% 1, 0, NA))

# remove everything but the cleaned data
rm(db0, dbP, Prices, ETH2013)

# 
db1 <- unique(db1)

# get square and interaction terms for translog
db1$logNsq <- db1$logN^2
db1$loglabsq <- db1$loglab^2
db1$logseedsq <- db1$logseed^2
db1$logarea_tot <- log(db1$area_tot)

# Also divide GGD and AI by 1000
# to get a more interpretable coef
db1$GGD <- db1$GGD/1000
db1$AI <- db1$AI/1000

# there are several core variables which
# have extreme values. 
db1 <- filter(db1,
              lab < 3000,
              seedha < 1000 | is.na(seedha)) 


saveRDS(db1, file.path(root, "Cache/db1.rds"))
