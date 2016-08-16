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

# Select maize fields per survey year
# Note that we only know on which field fertilizer is used, not if they are maize plots.
# We decide to calculate the average fertilizer price over maize plots only as it is possible that because of subsidies or other policies 
# the price of the same type of fertilizer (e.g. UREA) can differ between type of crop, even in the same region

ETH2011 <- filter(dbP, surveyyear ==2011) %>% select(-WPn, -crop_price)
ETH2013 <- filter(dbP, surveyyear ==2013) %>% select(-WPn, -crop_price)

# read in the fertilizer data, linkin location data and combine in one file

fert2011_1 <- remove_all_labels(read_dta(file.path(dataPath, "ETH/2011/Data/sect3_pp_w1.dta"))) %>%
  transmute(holder_id, household_id, parcel_id, field_id, typ=pp_s3q15,
            qty1=pp_s3q16_a, qty2=pp_s3q16_b/1000) %>%
  mutate(qty1 = ifelse(is.na(qty1), 0, qty1),
         qty2 = ifelse(is.na(qty2), 0, qty2),
         qty=qty1+qty2,
         typ = ifelse(typ %in% 1, "UREA", NA)) %>%
  select(-qty1, -qty2)

fert2011_1 <- remove_all_labels(read_dta(file.path(dataPath, "ETH/2011/Data/sect3_pp_w1.dta"))) %>%
  transmute(holder_id, household_id, parcel_id, field_id, typ=pp_s3q18,
            qty1=pp_s3q19_a, qty2=pp_s3q19_b/1000) %>%
  mutate(qty1 = ifelse(is.na(qty1), 0, qty1),
         qty2 = ifelse(is.na(qty2), 0, qty2),
         qty=qty1+qty2,
         typ = ifelse(typ %in% 1, "DAP", NA)) %>%
  select(-qty1, -qty2)

fert2011 <- bind_rows(fert2011_1, fert2013_1)


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

df <- fertmar
level <- "ZONE"
group <- "ZONE"

library(lazyeval)
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
# Values are winsored and aggregates are presented for at least 5 values (not a problem here)
maizemar <- filter(dbP, surveyyear ==2013) %>% 
  select(ZONE, REGNAME, WOREDACODE, KEBELECODE, price = crop_price) %>%
  mutate(price = winsor2(price))

mpCountry <- maizemar %>% 
  dplyr::summarize(price = median(price, na.rm=T)) %>%
  mutate(level = "country") 
mpCountry <- mutate(base, price = mpCountry$price,
                    level = "country")

mpZone <- medianPrice_f(maizemar, "zone", c("ZONE"))
mpRegion <- medianPrice_f(maizemar, "region", c("ZONE", "REGNAME"))
mpWoreda <- medianPrice_f(maizemar, "woreda", c("ZONE", "REGNAME", "WOREDACODE"))
mpKebele <- medianPrice_f(maizemar, "kebele", c("ZONE", "REGNAME", "WOREDACODE", "KEBELECODE"))

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


# combine price data
Prices <- rbind(fertPrice, maizePrice)
saveRDS(Prices, file = "Cache\\Prices_ETH.rds")
