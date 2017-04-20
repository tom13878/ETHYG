#'========================================================================================================================================
#' Project:  IMAGINE
#' Subject:  Script to create figures
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
p_load("frontier")

### DETERMINE ROOT PATH
root <- find_root(is_rstudio_project)

### DATAPATH

### R SETTINGS
options(scipen=999) # surpress scientific notation
options("stringsAsFactors"=FALSE) # ensures that characterdata that is loaded (e.g. csv) is not turned into factors
options(digits=4)

### SOURCE
source(file.path(root, "Code/waterfall_plot.r"))

### LOAD DATA
db3 <- readRDS(file.path(root, "Cache/db3.rds"))


# Group ZONES
db3 <- db3 %>% 
  mutate(ZONE = ifelse( ZONE %in% c("DIRE DAWA", "AFAR", "GAMBELLA"), "Other", ZONE))


# table of N, N0, pm, pn, relprice and mpp by zone
by_zone <- db3 %>%
  ungroup() %>%
  group_by(ZONE) %>%
  summarize(number = n(),
            n = mean(N, na.rm=TRUE),
            n0 = mean(N[yesN == 1], na.rm=TRUE),
            Npm = mean(Npm, na.rm=TRUE),
            pn = mean(Pn, na.rm=TRUE),
            pm = mean(Pm, na.rm=TRUE),
            relprice = mean(relprice, na.rm=TRUE),
            mpp = mean(as.numeric(mpp), na.rm=TRUE),
            vcr = mean(as.numeric(vcr), na.rm = TRUE))

# get rid of NA values for the yield gap. caused by NA values
# for the Npm variable
db4 <- filter(db3, !is.na(EY))


         
         
# Table with yield levels
yld_level <- bind_rows(
  db4 %>% 
    select(ZONE, Y, Ycor, TEY, EY, PFY, PY, area) %>%
    gather(yld_level, yld, -ZONE, -area) %>%
    group_by(yld_level, ZONE) %>%
    summarize(yld = (sum((yld)*area)/sum(area))) %>%
    spread(yld_level, yld),
  db4 %>% 
    select(ZONE, Y, Ycor, TEY, EY, PFY, PY, area) %>%
    gather(yld_level, yld, -ZONE, -area) %>%
    group_by(yld_level) %>%
    summarize(ZONE = "Total", yld = (sum((yld)*area)/sum(area))) %>%
    spread(yld_level, yld)) %>%
  select(ZONE, Y, Ycor, TEY, EY, PFY, PY)


# Table with absolute yield gap information per zone
# Note that by definition, YG_s computed by weighting individual YG_s values is not the same as multiplication of weighted TEYG_s etc.
# We therefore calculate YG_s as the product of the weighted components.
yg_level <- bind_rows(
  db4 %>% 
    dplyr::select(ZONE, TEYG_l, EYG_l, EUYG_l, TYG_l, YG_l_Ycor, YG_l, area) %>%
    gather(yg_level, yg, -ZONE, -area) %>%
    group_by(yg_level, ZONE) %>%
    summarize(yg = (sum((yg)*area)/sum(area))) %>%
    spread(yg_level, yg),
  db4 %>% 
    dplyr::select(ZONE, TEYG_l, EYG_l, EUYG_l, TYG_l, YG_l_Ycor, YG_l, area) %>%
    gather(yg_level, yg, -ZONE, -area) %>%
    group_by(yg_level) %>%
    summarize(ZONE = "Total", yg = (sum((yg)*area)/sum(area))) %>%
    spread(yg_level, yg))
  

# Table with yield gap shares
yg_share <- yg_level %>%
  mutate(
    TEYG = 100*TEYG_l/YG_l_Ycor,
    EYG = 100*EYG_l/YG_l_Ycor,
    EUYG = 100*EUYG_l/YG_l_Ycor,
    TYG = 100*TYG_l/YG_l_Ycor,
    YG = 100*(TEYG_l + EYG_l + EUYG_l + TYG_l)/YG_l_Ycor) %>%
  dplyr::select(ZONE, TEYG:YG) 

