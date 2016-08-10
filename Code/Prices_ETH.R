#######################################
###### ANALYSIS of ETH price data #####
#######################################

# CHECK
# Median of prices
# Compute community prices!
# Compare prices with other prices and check if they are realistic!

wdPath <- "D:\\Data\\Projects\\ETHYG"
extraDataPath <- "D:\\Dropbox\\Michiel_research\\2285000066 Africa Maize Yield Gap\\Analysis\\ETH\\Data"
setwd(wdPath)

surveyPath <- "N:\\Internationaal Beleid  (IB)\\Projecten\\2285000066 Africa Maize Yield Gap\\SurveyData"

library(plyr)
library(dplyr)
library(ggplot2)
library(stargazer)
library(haven)
library(tidyr)
library(xtable)
library(DescTools)
library(sandwich)
library(lmtest)
library(assertive)

options(scipen=999)

# winsor code
source("Code/winsor.R")

# read in the base data and select relevant variables
ETH_2013 <- readRDS("data/ETH_data_2013.rds") %>%
              select(holder_id:field_id, ea_id2, lon, lat, region, district, region_lsms, AEZ)


# read in the fertilizer data, linkin location data and combine in one file
fert2013_1 <- read_dta(file.path(surveyPath, "ETH\\2013\\Data\\Post-Planting\\sect3_pp_w2.dta")) %>%
              dplyr::select(holder_id, household_id2, parcel_id, field_id, typ=pp_s3q15, qty=pp_s3q16_a,
                purch=pp_s3q16b, purch_kg=pp_s3q16c, valu=pp_s3q16d) %>%
  filter(!(household_id2 %in% "")) %>%
  mutate(typ = ifelse(typ %in% 1, "UREA", NA),
         purch = ifelse(purch %in% 1, 1, 0))

# WDswitch  
fert2013_2 <- read_dta(file.path(surveyPath, "ETH\\2013\\Data\\Post-Planting\\sect3_pp_w2.dta")) %>%
  dplyr::select(holder_id, household_id2, parcel_id, field_id, typ=pp_s3q18, qty=pp_s3q19_a,
                purch=pp_s3q19b, purch_kg=pp_s3q19c, valu=pp_s3q19d) %>%
  filter(!(household_id2 %in% "")) %>%
  mutate(typ = ifelse(typ %in% 1, "DAP", NA),
       purch = ifelse(purch %in% 1, 1, 0))

# Combine fertilizer files and join with RTH_2013 to select only maize plots. 
fert2013 <- rbind(fert2013_1, fert2013_2)
fert2013[] <- lapply(fert2013, strip_attributes)
fert2013 <- fert2013 %>% 
  left_join(ETH_2013,.) %>%
  unique() %>% 
  do(filter(., complete.cases(.)))
  
  

# -------------------------------------
# read in nitrogen conversion file

conv <- read.csv(file.path(extraDataPath, "Fert_comp.csv")) %>%
  transmute(typ=Fert_type2, n=N_share/100, p=P_share/100) %>%
  filter(typ %in% c("UREA", "DAP"))

fert <- left_join(fert2013, conv)

fert <- mutate(fert,
               Vfert=valu/purch_kg,
               Qn=qty*n,
               Qp=qty*p,
               price=Vfert/n)

# Distribution of fertilizer types
# Find out which combinations are used and add as potential regressor.

# Recode AEZ into 3 zones
fert$AEZ2 <- fert$AEZ
fert$AEZ2 <- mapvalues(fert$AEZ2, from = c("Tropic-warm / semi-arid"), to = c("Tropic-warm"))
fert$AEZ2 <- mapvalues(fert$AEZ2, from = c("Tropic-warm / sub-humid"), to = c("Tropic-warm"))
fert$AEZ2 <- factor(fert$AEZ2)

fertdist <- fert %>% group_by(AEZ2, typ) %>%
  filter(!is.na(typ)) %>%
  summarise(number = n())


# Construct price data.frame
# Note that some ea_id2 codes have two or more lat, lon coordinates, which should not be possible. For this reason I group by lat, lon, not ea_id2
# construct base dataframe with all zones and regions
base <- ETH_2013 %>% 
  filter(!is.na(region_lsms)) %>%
  dplyr::select(region_lsms, district, lat, lon) %>%
  unique()
  
# market  prices
fertmar <- fert %>%
  mutate(price = winsor2(price))

pricesCountry <- fertmar %>% 
  dplyr::summarize(
    number = sum(!is.na(price)),
    price = mean(price, na.rm=T)) %>%
  mutate(level = "country")

pricesPerZone <- fertmar %>% # note that I use region as zone because regions are relatively large in ETH. I use district as region.
  group_by(region_lsms) %>%
  dplyr::summarize(
    number = sum(!is.na(price)),
    price = mean(price, na.rm=T)) %>%
  filter(number>=5) %>%
  mutate(level = "zone")

pricesPerRegion <- fertmar %>% 
  group_by(district) %>%
  dplyr::summarize(
    number = sum(!is.na(price)),
    price = mean(price, na.rm=T)) %>%
  filter(number>=5) %>%
  mutate(level = "district")

pricesPerEA <- fertmar %>% 
  group_by(lon, lat) %>%
  dplyr::summarize(
    number = sum(!is.na(price)),
    price = mean(price, na.rm=T)) %>%
  filter(number>=5) %>%
  mutate(level = "ea")

pricesRegionInter <- fertmar %>%
  group_by(lon, lat, district, region_lsms) %>%
  dplyr::summarize(
    number = sum(!is.na(price)),
    price = mean(price, na.rm=T)) %>%
  left_join(base, .) %>%
  filter(is.na(price) | number<5) %>%
  dplyr::select(-number, -price) %>%
  left_join(., pricesPerRegion) %>%
  filter(!is.na(price))

pricesZoneInter <- fertmar %>%
  group_by(lon, lat, district, region_lsms) %>%
  dplyr::summarize(
    number = sum(!is.na(price)),
    price = mean(price, na.rm=T)) %>%
  left_join(base, .) %>%
  filter(is.na(price) | number<5) %>%
  dplyr::select(-number, -price) %>%
  left_join(., pricesPerRegion) %>%
  filter(is.na(price)) %>%
  dplyr::select(-number, -price, -level) %>%
  left_join(., pricesPerZone) %>%
  filter(!is.na(price))

pricesCountryInter <- fertmar %>%
  group_by(lon, lat, district, region_lsms) %>%
  dplyr::summarize(
    number = sum(!is.na(price)),
    price = mean(price, na.rm=T)) %>%
  left_join(base, .) %>%
  filter(is.na(price) | number<5) %>%
  dplyr::select(-number, -price) %>%
  left_join(., pricesPerRegion) %>%
  filter(is.na(price)) %>%
  dplyr::select(-number, -price, -level) %>%
  left_join(., pricesPerZone) %>%
  filter(is.na(price)) %>%
  dplyr::select(-number, -price, -level) %>%
  mutate(number = pricesCountry$number,
         price = pricesCountry$price,
         level = pricesCountry$level)

pricesPerEAmar <- fertmar %>%
  group_by(lon, lat, district, region_lsms) %>%
  dplyr::summarize(
    number = sum(!is.na(price)),
    price = mean(price, na.rm=T)) %>%
  mutate(level = "ea_id") %>%
  left_join(base, .) %>%
  filter(!is.na(price) & number>=5) %>%
  rbind(., pricesCountryInter, pricesZoneInter,pricesRegionInter) %>%
  ungroup() %>%
  mutate(product = "fertilizer", type = "Pn")

rm(pricesCountry, pricesPerZone, pricesPerRegion, pricesPerEA, pricesCountryInter, pricesRegionInter,pricesZoneInter)

# Maize prices
# Values are winsored and aggregates are presented for at least 5 values (not a problem here)
ETH_2013 <- select(ETH_2013, holder_id, household_id2, region_lsms, lat, lon, ea_id2, region, district) %>%
          unique()
maizePrices <- read_dta(file.path(surveyPath, "ETH\\2013\\Data/Post-Harvest/sect11_ph_w2.dta")) %>%
  select(holder_id, household_id2, crop_name, sold=ph_s11q01, qty_sold_kg=ph_s11q03_a,
         valu =ph_s11q04) %>%
  filter(crop_name %in% "MAIZE" & sold==1) %>%
  mutate(price = valu/qty_sold_kg)
  
maizePrices[] <- lapply(maizePrices, strip_attributes)

maizePrices <- left_join(maizePrices, ETH_2013)

mpricesCountry <- maizePrices %>% 
  dplyr::summarize(
    number = sum(!is.na(price)),
    price = mean(price, na.rm=T)) %>%
  mutate(level = "country")

mpricesPerZone <- maizePrices %>% # note that I use region as zone because regions are relatively large in ETH. I use district as region.
  group_by(region_lsms) %>%
  dplyr::summarize(
    number = sum(!is.na(price)),
    price = mean(price, na.rm=T)) %>%
  filter(number>=5) %>%
  mutate(level = "zone")

mpricesPerRegion <- maizePrices %>% 
  group_by(district) %>%
  dplyr::summarize(
    number = sum(!is.na(price)),
    price = mean(price, na.rm=T)) %>%
  filter(number>=5) %>%
  mutate(level = "district")

mpricesPerEA <- maizePrices %>% 
  group_by(lon, lat) %>%
  dplyr::summarize(
    number = sum(!is.na(price)),
    price = mean(price, na.rm=T)) %>%
  filter(number>=5) %>%
  mutate(level = "ea")

mpricesRegionInter <- maizePrices %>%
  group_by(lon, lat, district, region_lsms) %>%
  dplyr::summarize(
    number = sum(!is.na(price)),
    price = mean(price, na.rm=T)) %>%
  left_join(base, .) %>%
  filter(is.na(price) | number<5) %>%
  dplyr::select(-number, -price) %>%
  left_join(., mpricesPerRegion) %>%
  filter(!is.na(price))

mpricesZoneInter <- maizePrices %>%
  group_by(lon, lat, district, region_lsms) %>%
  dplyr::summarize(
    number = sum(!is.na(price)),
    price = mean(price, na.rm=T)) %>%
  left_join(base, .) %>%
  filter(is.na(price) | number<5) %>%
  dplyr::select(-number, -price) %>%
  left_join(., mpricesPerRegion) %>%
  filter(is.na(price)) %>%
  dplyr::select(-number, -price, -level) %>%
  left_join(., mpricesPerZone) %>%
  filter(!is.na(price))

mpricesCountryInter <- maizePrices %>%
  group_by(lon, lat, district, region_lsms) %>%
  dplyr::summarize(
    number = sum(!is.na(price)),
    price = mean(price, na.rm=T)) %>%
  left_join(base, .) %>%
  filter(is.na(price) | number<5) %>%
  dplyr::select(-number, -price) %>%
  left_join(., mpricesPerRegion) %>%
  filter(is.na(price)) %>%
  dplyr::select(-number, -price, -level) %>%
  left_join(., mpricesPerZone) %>%
  filter(is.na(price)) %>%
  dplyr::select(-number, -price, -level) %>%
  mutate(number = mpricesCountry$number,
         price = mpricesCountry$price,
         level = mpricesCountry$level)

maizePricesPerEA <- maizePrices %>%
  group_by(lon, lat, district, region_lsms) %>%
  dplyr::summarize(
    number = sum(!is.na(price)),
    price = mean(price, na.rm=T)) %>%
  mutate(level = "ea_id") %>%
  left_join(base, .) %>%
  filter(!is.na(price) & number>=5) %>%
  rbind(., mpricesCountryInter, mpricesZoneInter, mpricesRegionInter) %>%
  ungroup() %>%
  mutate(product = "maize", type = "Pm")

rm(mpricesCountry, mpricesPerZone, mpricesPerRegion, mpricesPerEA, mpricesRegionInter, mpricesZoneInter, mpricesCountryInter)

# combine price data
Prices <- rbind(pricesPerEAmar, maizePricesPerEA) %>%
          mutate(region = as.factor(district)) %>%
          droplevels(.)

saveRDS(Prices, file = "Data\\Prices_ETH.rds")
