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

# Load pooled data
dbP <- readRDS(file.path(root, "Cache/Pooled_ETH.rds"))

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
    N = N/area)


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
  dplyr::select(hhid, ea_id, DISNAME = ZONENAME, REGNAME, parcel_id, field_id, holder_id,
                AEZ, fs,
                rain_year, rain_wq, SPEI,
                YA, YW, YP,
                slope, elevation,
                nutr_av,
                yld, 
                harv_lab, harv_lab_hire ,
                impr, 
                fung, herb,
                N, P, 
                manure, compost, other_org,
                crop_stand, cropping,
                legume, irrig, 
                area, area_tot, area_gps,
                sex, age,
                literate, cage, ed_any, N1555, family_size, death,
                dist_hh, dist_road, dist_market, dist_popcenter, dist_regcap,
                title,
                popEA,
                road, cost2small_town, bank, micro_finance, ext_agent, extension, 
                crop_count, surveyyear,
                rural, 
                lat, lon)


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
                       logarea = log(area_gps), # area_gps not area because we want to add plot size as proxy for economies of scale
                       rain_wq2 = rain_wq*rain_wq,
                       rain_year2 = rain_year*rain_year,
                       lograin = log(rain_year),
                       sex = as.numeric(ifelse(sex == "MALE", 0, 1)),
                       surveyyear2 = replace(surveyyear==2011, 1, 0))

db0 <- droplevels(db0)

# Load price data 
Prices <- readRDS(file.path(root, "cache/Prices_ETH.rds")) %>% select(hhid, fertilizer, maize)

# Merge with panel data
db1 <- left_join(db0, Prices) %>%
  mutate(relprice = fertilizer/maize) %>%
  rename(Pm = maize, Pn = fertilizer)

# Drop unused levels (e.g. Zanzibar in zone), which are giving problems with sfa
db1 <- droplevels(db1)

# the crop stand variable has problems
db1$crop_stand <- ifelse(db1$crop_stand %in% c("NAN", "OTHER LAND USE (SPECIFY)"), NA, db1$crop_stand)
db1$crop_stand <- ifelse(db1$crop_stand %in% c("PURESTAND"), "PURE STAND", db1$crop_stand)

# remove everything but the cleaned data
rm(db0, dbP, Prices)
