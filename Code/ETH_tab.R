#'========================================================================================================================================
#' Project:  IMAGINE
#' Subject:  Script to create tables
#' Author:   Michiel van Dijk
#' Contact:  michiel.vandijk@wur.nl
#'========================================================================================================================================

### PACKAGES
if(!require(pacman)) install.packages("pacman")
# Key packages
p_load("tidyverse", "readxl", "stringr", "scales", "RColorBrewer", "rprojroot")
# Spatial packages
#p_load("rgdal", "ggmap", "raster", "rasterVis", "rgeos", "sp", "mapproj", "maptools", "proj4", "gdalUtils")
# Additional packages
p_load("WDI", "countrycode", "stargazer")

### DETERMINE ROOT PATH
root <- find_root(is_rstudio_project)

### DATAPATH
source(file.path(root, "Code/get_dataPath_ETH.R"))
### R SETTINGS
options(scipen=999) # surpress scientific notation
options("stringsAsFactors"=FALSE) # ensures that characterdata that is loaded (e.g. csv) is not turned into factors
options(digits=4)

### SOURCE
source(file.path(root, "Code/sfaTable.R"))

### LOAD DATA
db9 <- readRDS(file.path(root, "Cache/db9.rds"))
#db_sfaCD_CRE_Z <- readRDS(file.path(root, "Cache/db_sfaCD_CRE_Z.rds"))
db1 <- readRDS(file.path(root, "Cache/db1.rds"))

# SUMMARY STATISTICS
# Number of unique households
length(unique(db9$hhid))

# summary statistics
# dbsum <- db_sfaCD_CRE_Z %>% dplyr::select(Yield = yld, ImprovedSeeds = impr, Slope = slope, 
#                                           yesN, Nitrogen = N, Rain = rain_wq, Irrigation = irrig, SOC = SOC2, phdum, 
#                                           Labour = lab,  Area = area, crop_count2, AEZ, surveyyear2) 

# Stargazer does not work with rmarkdown when pushed to word
# Need to create a function or code that creates a dataframe with summary statistics.
#stargazer(as.data.frame(dbsum), type = "text", digits=2, out=".\\FigTab\\summaryStat.html")

# Table with key information per zone and subtotals
Yieldsum <- bind_rows(
  db9 %>% 
    dplyr::select(Y, N, yesN, ZONE, surveyyear, area) %>%
    group_by(ZONE) %>%
    summarize(Yield = mean(Y),
              Yield_w = (sum(Y*area)/sum(area)),
              NitrogenUser = round(mean(yesN)*100, digits=1),
              Number=n()),
  db9 %>%
    dplyr::select(Y, N, yesN, area, ZONE) %>%
    summarize(ZONE = "Total",  
              Yield = mean(Y),
              Yield_w = (sum(Y*area)/sum(area)),
              NitrogenUser = round(mean(yesN)*100, digits=1),
              Number=n())
) %>% arrange(ZONE)

# Table with key information per zone and subtotals for pure maize plots
Yieldsum_pure <- bind_rows(
  db9 %>%
    filter(crop_count2==1) %>%
    dplyr::select(Y, N, yesN, ZONE, surveyyear, area) %>%
    group_by(ZONE) %>%
    summarize(Yield_p = mean(Y),
              Yield_w_p = (sum(Y*area)/sum(area))
    ),
  db9 %>%
    filter(crop_count2==1) %>%
    dplyr::select(Y, N, yesN, ZONE, area) %>%
    summarize(ZONE = "Total", 
              Yield_p = mean(Y),
              Yield_w_p = (sum(Y*area)/sum(area))
    )
)



Nitrogensum <- bind_rows(
  db9 %>% 
    dplyr::select(Y, N, yesN, ZONE) %>%
    filter(yesN ==1) %>%
    group_by(ZONE) %>%
    summarize(Nitrogen = mean(N)),
  db9 %>%
    dplyr::select(Y, N, yesN, ZONE) %>%
    filter(yesN ==1) %>%
    summarize(ZONE= "Total", Nitrogen = mean(N))
) 

Pricessum <- bind_rows(
  db1 %>% 
    dplyr::select(ZONE, Pn, Pm) %>% 
    do(filter(., complete.cases(.))) %>%
    group_by(ZONE) %>%
    summarize(
      NitrogenPrice = mean(Pn, na.rm=T),
      MaizePrice = round(mean(Pm, na.rm=T), digits=0)),
  db1 %>% 
    dplyr::select(ZONE, Pn, Pm) %>% 
    do(filter(., complete.cases(.))) %>%
    summarize(ZONE = "Total",
              NitrogenPrice = mean(Pn, na.rm=T),
              MaizePrice = round(mean(Pm, na.rm=T), digits=0))
) 


Zonalsum <- left_join(Yieldsum, Nitrogensum) %>%
  left_join(., Yieldsum_pure) %>%
  left_join(., Pricessum) %>%
  dplyr::select(ZONE, Number, Yield_w, Yield_w_p, NitrogenUser, Nitrogen, NitrogenPrice, MaizePrice) %>%
  mutate(ZONE = factor(ZONE, levels = c("TIGRAY", "AMHARA", "OROMIYA", "SOMALI", "BENSHANGULGUMUZ", "SNNP", "HARARI", "GAMBELLA", "DIRE DAWA", "Total"))) %>%
  arrange(ZONE)


# SFA table
# sfaCD_CRE_Z <- sfa(logyld ~ noN + logN + loglab + dumoxen + 
#                      logarea + irrig + impr + slope + elevation + SOC2 + phdum2 + 
#                      rain_wq + AEZ +
#                      crop_count2 + surveyyear2 + 
#                      noN_bar + logN_bar + loglab_bar + logarea_bar + oxen_bar +
#                      irrig_bar + impr_bar + slope_bar + rain_wq_bar + crop_count_bar + 
#                      r | 
#                      sex + age +  title + literate + ed_any + extension + credit +
#                      dist_hh + dist_market + popEA -1
#                    ,data = db1, maxit = 1500, restartMax = 20, tol = 0.000001)
# summary(sfaCD_CRE_Z, extraPar = TRUE)
# lrtest(sfaCD_CRE_Z)
# 
# xtable <- sfaTable_f(sfaCD_CRE_Z)[[1]]
# ztable <- sfaTable_f(sfaCD_CRE_Z)[[2]]


# Table with yield levels
YieldLevels <- bind_rows(
  db9 %>% 
    dplyr::select(Zone = ZONE, Y, Ycor, TEY, EY, PFY, PY, area) %>%
    group_by(Zone) %>%
    summarize(Y =(sum((Y)*area)/sum(area)),
              Ycor = (sum((Ycor)*area)/sum(area)),
              TEY = (sum((TEY)*area)/sum(area)),
              EY = (sum((EY)*area, na.rm=TRUE)/sum(area)),
              PFY = (sum((PFY)*area)/sum(area)),
              PY = (sum((PY)*area)/sum(area))
    ),
  db9 %>% 
    dplyr::select(Zone = ZONE, Y, Ycor, TEY, EY, PFY, PY, area) %>%
    summarize(Zone = "Total", 
              Y =(sum((Y)*area)/sum(area)),
              Ycor = (sum((Ycor)*area)/sum(area)),
              TEY = (sum((TEY)*area)/sum(area)),
              EY = (sum((EY)*area, na.rm=TRUE)/sum(area)),
              PFY = (sum((PFY)*area)/sum(area)),
              PY = (sum((PY)*area)/sum(area)))) %>%
  dplyr::select(Zone, Y, Ycor, TEY, EY, PFY, PY)


# Table with relative yield gap information per zone
# Note that by definition, YG_s computed by weighting individual YG_s values is not the same as multiplication of weighted TEYG_s etc.
# We therefore calculate YG_s as the product of the weighted components.
ZonalYieldGap_s <- bind_rows(
  db9 %>% 
    dplyr::select(Zone = ZONE, ERROR_s, TEYG_s, EYG_s, EUYG_s, TYG_s, YG_s_Ycor, YG_s, area) %>%
    group_by(Zone) %>%
    summarize(ERROR_s =(sum((ERROR_s)*area, na.rm=TRUE)/sum(area)),
              TEYG_s = (sum((TEYG_s)*area, na.rm=TRUE)/sum(area)),
              EYG_s = (sum((EYG_s)*area, na.rm=TRUE)/sum(area)),
              EUYG_s = (sum((EUYG_s)*area, na.rm=TRUE)/sum(area)),
              TYG_s = (sum((TYG_s)*area, na.rm=TRUE)/sum(area)),
              YG_s = (sum((YG_s)*area, na.rm=TRUE)/sum(area)),
              YG_s_Ycor = (sum((YG_s_Ycor)*area, na.rm=TRUE)/sum(area)),
              ERROR = (1-ERROR_s)*100,
              TEYG = (1-TEYG_s)*100,
              EYG = (1-EYG_s)*100,
              EUYG = (1-EUYG_s)*100,
              TYG = (1-TYG_s)*100,
              YG = (1-(ERROR_s*TEYG_s*EYG_s*EUYG_s*TYG_s))*100,
              YG_Ycor = (1-(TEYG_s*EYG_s*EUYG_s*TYG_s))*100
    ),
  db9 %>% 
    dplyr::select(Zone = ZONE, ERROR_s, TEYG_s, EYG_s, EUYG_s, TYG_s, YG_s_Ycor, YG_s, area) %>%
    summarize(Zone = "Total", 
              ERROR_s =(sum((ERROR_s)*area, na.rm=TRUE)/sum(area)),
              TEYG_s = (sum((TEYG_s)*area, na.rm=TRUE)/sum(area)),
              EYG_s = (sum((EYG_s)*area, na.rm=TRUE)/sum(area)),
              EUYG_s = (sum((EUYG_s)*area, na.rm=TRUE)/sum(area)),
              TYG_s = (sum((TYG_s)*area, na.rm=TRUE)/sum(area)),
              YG_s = (sum((YG_s)*area, na.rm=TRUE)/sum(area)),
              YG_s_Ycor = (sum((YG_s_Ycor)*area, na.rm=TRUE)/sum(area)),
              ERROR = (1-ERROR_s)*100,
              TEYG = (1-TEYG_s)*100,
              EYG = (1-EYG_s)*100,
              EUYG = (1-EUYG_s)*100,
              TYG = (1-TYG_s)*100,
              YG = (1-(ERROR_s*TEYG_s*EYG_s*EUYG_s*TYG_s))*100,
              YG_Ycor = (1-(TEYG_s*EYG_s*EUYG_s*TYG_s))*100)) %>%
  dplyr::select(Zone, TEYG, EYG, EUYG, TYG, YG = YG_Ycor)


