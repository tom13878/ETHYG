#######################################
###### ANALYSIS of ETH price data #####
#######################################

# Compare prices with other prices and check if they are realistic!

wdPath <- "D:\\Data\\Projects\\ETHYG"
setwd(wdPath)

dataPath <- "N:\\Internationaal Beleid  (IB)\\Projecten\\2285000066 Africa Maize Yield Gap\\SurveyData"

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
library(sjmisc)
library(lazyeval)

options(scipen=999)

# winsor code
source("Code/winsor.R")

#######################################
############## LOAD DATA ##############
#######################################

# Load pooled data 
dbP <- readRDS("Cache/Pooled_ETH.rds") %>%
 filter(crop_code %in% 2 & status %in% "HEAD") %>%
 select(hhid, ZONE = REGNAME, REGNAME = ZONENAME, WOREDACODE, KEBELECODE, parcel_id, field_id, holder_id, # ZONE AND REGNAMES reversed to remain consistent with other LSMS
        WPn, crop_price, surveyyear) 

# read in nitrogen conversion file
conv <- read.csv(file.path(dataPath, "Other/Fertilizer/Fert_comp.csv")) %>%
  transmute(typ=Fert_type2, n=N_share/100, p=P_share/100) %>%
  filter(typ %in% c("UREA", "DAP"))

# Note that we only know on which field fertilizer is used, not if they are maize plots.
# We decide to calculate the average fertilizer price over maize plots only as it is possible that because of subsidies or other policies 
# the price of the same type of fertilizer (e.g. UREA) can differ between type of crop, even in the same region
# There is no information on kg and value of fertilizer purchased for 2011. We will use 2013 data for 2011.
# There is information  on maize prices for 2011, which we will use at plot level. However due to the small number of observations
# we will use the regional averages from 2013 but correct for inflation.  
# Inflation is not a problem if we use relative prices (pmaize/pfert), which are assumed to have same inflation.

# Select maize fields per survey year
ETH2013 <- filter(dbP, surveyyear ==2013) %>% 
  select(-WPn, -crop_price) %>% 
  unique # There are four hhid-parcel combinations that have two entries, one with price data, one without => duplicates deleted

# read in the fertilizer data, linkin location data and combine in one file
fert2013_1 <- read_dta(file.path(dataPath, "ETH\\2013\\Data\\Post-Planting\\sect3_pp_w2.dta")) %>%
  dplyr::select(holder_id, household_id, household_id2, parcel_id, field_id, typ=pp_s3q15, purch_kg=pp_s3q16c, valu=pp_s3q16d) %>%
  mutate(household_id = zap_empty(household_id),
         typ = ifelse(typ %in% 1, "UREA", NA))

fert2013_2 <- read_dta(file.path(dataPath, "ETH\\2013\\Data\\Post-Planting\\sect3_pp_w2.dta")) %>%
  dplyr::select(holder_id, household_id, household_id2, parcel_id, field_id, typ=pp_s3q18, purch_kg=pp_s3q19c, valu=pp_s3q19d) %>%
  mutate(household_id = zap_empty(household_id),
         typ = ifelse(typ %in% 1, "DAP", NA))

# Combine fertilizer files and join with ETH2013 to select only maize plots. 
fert2013 <- remove_all_labels(rbind(fert2013_1, fert2013_2)) %>%
            mutate(hhid = ifelse(is.na(household_id), household_id2, household_id)) %>% # Ensure that hhid is the same as in dbP
            dplyr::select(-household_id, -household_id2) %>%
            do(filter(., complete.cases(.))) %>%
            left_join(ETH2013,.) %>%
            unique() %>% 
            do(filter(., complete.cases(.))) %>%
            left_join(., conv) %>%
            mutate(purch_kg = ifelse(purch_kg == 0, NA, purch_kg),
                   valu = ifelse(valu == 0, NA, valu),
                   Qn=purch_kg*n,
                   price = valu/Qn) 

# Construct price data.frame
# Note that some ea_id2 codes have two or more lat, lon coordinates, which should not be possible. For this reason we use KEBELE as lowest level for calculating median prices
# construct base dataframe with all zones, regions, woreda and kebele

base <- dbP %>% 
  dplyr::select(ZONE, REGNAME, WOREDACODE, KEBELECODE) %>%
  unique() %>%
  na.omit

# Values are winsored aggregates are presented for at least 5 values
# market  prices
fertmar <- fert2013 %>%
  mutate(price = winsor2(price))

medianPrice_f <- function(df, level, group){
  prices <- df %>% 
    group_by_(.dots = c(group)) %>%
    dplyr::summarize(
      number = sum(!is.na(price)),
      price = median(price, na.rm=T)) %>%
    filter(number>=5) %>%
    mutate(level = level) %>%
    select(-number) 
   #prices <- setNames(prices, c(group, "price", "level")) 
   out <- left_join(base, prices)
   return(out)
}

fpCountry <- fertmar %>% 
  dplyr::summarize(price = median(price, na.rm=T)) %>%
  mutate(level = "country") 
fpCountry <- mutate(base, price = fpCountry$price,
                       level = "country")

fpZone <- medianPrice_f(fertmar, "zone", c("ZONE"))
fpRegion <- medianPrice_f(fertmar, "region", c("ZONE", "REGNAME"))
fpWoreda <- medianPrice_f(fertmar, "woreda", c("ZONE", "REGNAME", "WOREDACODE"))
fpKebele <- medianPrice_f(fertmar, "kebele", c("ZONE", "REGNAME", "WOREDACODE", "KEBELECODE"))

fertPrice <- bind_rows(fpWoreda, fpKebele, fpRegion, fpZone, fpCountry) %>%
  na.omit %>%
  spread(level, price) %>%
  mutate(price = ifelse(!is.na(kebele), kebele, 
                  ifelse(!is.na(woreda), woreda,
                         ifelse(!is.na(zone), zone, country))),
         source = ifelse(!is.na(kebele), "kebele", 
                        ifelse(!is.na(woreda), "woreda",
                               ifelse(!is.na(zone), "zone", "country"))),
         product = "fertilizer") %>%
  select(-country, -zone, -region, -woreda, -kebele)

# Maize prices
# Using inflation rate for 2011 and 2012. These years were selected as the main part of the survey takes place in these years.
# from world bank:
# http://data.worldbank.org/indicator/FP.CPI.TOTL.ZG/countries/TZ?display=graph
# -------------------------------------

inflation <- read.csv(file.path(dataPath,"Other/Inflation/inflation.csv"))
rate2012 <- inflation$inflation[inflation$code=="ET" & inflation$year==2012]/100
rate2013 <- inflation$inflation[inflation$code=="ET" & inflation$year==2013]/100
inflate <- (1 + rate2012)*(1 + rate2013)

# Values are winsored per surveyyear. 2011 values are inflated to 2013 levels and used at hhid level. 
# Regional values for 2013 are used to impute missing values for 2011 and 2013.
# Aggregates are presented for at least 5 values.
maizemar2011 <- filter(dbP, surveyyear ==2011) %>% 
  select(hhid, holder_id, field_id, parcel_id, ZONE, REGNAME, WOREDACODE, KEBELECODE, surveyyear, price = crop_price) %>%
  mutate(price = winsor2(price),
         price = price*inflate)

maizemar2013 <- filter(dbP, surveyyear ==2013) %>% 
  select(hhid, holder_id, field_id, parcel_id, ZONE, REGNAME, WOREDACODE, KEBELECODE, surveyyear, price = crop_price) %>%
  mutate(price = winsor2(price))

mpCountry <- maizemar2013 %>% 
  dplyr::summarize(price = median(price, na.rm=T)) %>%
  mutate(level = "country") 
mpCountry <- mutate(base, price = mpCountry$price,
                    level = "country")

mpZone <- medianPrice_f(maizemar2013, "zone", c("ZONE"))
mpRegion <- medianPrice_f(maizemar2013, "region", c("ZONE", "REGNAME"))
mpWoreda <- medianPrice_f(maizemar2013, "woreda", c("ZONE", "REGNAME", "WOREDACODE"))
mpKebele <- medianPrice_f(maizemar2013, "kebele", c("ZONE", "REGNAME", "WOREDACODE", "KEBELECODE"))

maizePrice <- bind_rows(mpWoreda, mpKebele, mpRegion, mpZone, mpCountry) %>%
  na.omit %>%
  spread(level, price) %>%
  mutate(price = ifelse(!is.na(kebele), kebele, 
                        ifelse(!is.na(woreda), woreda,
                               ifelse(!is.na(zone), zone, country))),
         source = ifelse(!is.na(kebele), "kebele", 
                         ifelse(!is.na(woreda), "woreda",
                                ifelse(!is.na(zone), "zone", "country"))),
         product = "maize") %>%
  select(-country, -zone, -region, -woreda, -kebele)

# Create price file at plot level and substitute regional prices when plot level price is not available.
# Again, we winsor the prices first
plotFertPrice <- select(dbP, hhid, holder_id, field_id, parcel_id, ZONE, REGNAME, WOREDACODE, KEBELECODE, surveyyear, WPn) %>%
    left_join(fertPrice) %>%
    mutate(WPn = winsor2(WPn),
           price = ifelse(is.na(WPn), price, WPn)) %>%
  select(-WPn) 
 
plotMaizePrice <- bind_rows(maizemar2011, maizemar2013) %>%
  rename(crop_price = price) %>%
  left_join(maizePrice) %>%
  mutate(price = ifelse(is.na(crop_price), price, crop_price)) %>%
  select(-crop_price)  

Prices <- rbind(plotFertPrice, plotMaizePrice) %>% 
  unique() %>% # There are four hhid-parcel combinations that have two entries, one with price data, one without => duplicates deleted
  select(-source) %>%
  spread(product, price)

saveRDS(Prices, file = "Cache\\Prices_ETH.rds")
