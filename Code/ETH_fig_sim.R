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
dbsim <- readRDS(file.path(root, "Cache/dbsim.rds"))

# get rid of NA values for the yield gap. caused by NA values
# for the Npm variable
db4 <- filter(dbsim, !is.na(EY))

# Table with yield levels
yld_level <- bind_rows(
  db4 %>% 
    select(ZONE, Y, Ycor, Yext, TEY, EY, Yseed, Yimpr, PFY, PY, area) %>%
    gather(yld_level, yld, -ZONE, -area) %>%
    group_by(yld_level, ZONE) %>%
    summarize(yld = (sum((yld)*area)/sum(area))) %>%
    spread(yld_level, yld),
  db4 %>% 
    select(ZONE, Y, Ycor, Yext, TEY, EY, Yseed, Yimpr, PFY, PY, area) %>%
    gather(yld_level, yld, -ZONE, -area) %>%
    group_by(yld_level) %>%
    summarize(ZONE = "Total", yld = (sum((yld)*area)/sum(area))) %>%
    spread(yld_level, yld)) %>%
  select(ZONE, Y, Ycor, Yext, TEY, EY, Yseed, Yimpr, PFY, PY)


# Table with absolute yield gap information per zone
# Note that by definition, YG_s computed by weighting individual YG_s values is not the same as multiplication of weighted TEYG_s etc.
# We therefore calculate YG_s as the product of the weighted components.
yg_level <- bind_rows(
  db4 %>% 
    dplyr::select(ZONE, TEYG_l, EYG_l, EUYG_l, TYG_l, YG_l_Ycor, YG_l, YextG_l, YseedG_l, YimprG_l, area) %>%
    gather(yg_level, yg, -ZONE, -area) %>%
    group_by(yg_level, ZONE) %>%
    summarize(yg = (sum((yg)*area)/sum(area))) %>%
    spread(yg_level, yg),
  db4 %>% 
    dplyr::select(ZONE, TEYG_l, EYG_l, EUYG_l, TYG_l, YG_l_Ycor, YG_l, YextG_l, YseedG_l, YimprG_l, area) %>%
    gather(yg_level, yg, -ZONE, -area) %>%
    group_by(yg_level) %>%
    summarize(ZONE = "Total", yg = (sum((yg)*area)/sum(area))) %>%
    spread(yg_level, yg))
  

# Table with yield gap shares
yg_share <- yg_level %>%
  mutate(
    TEYG = 100*TEYG_l/YG_l_Ycor,
    YextG = 100*YextG_l/YG_l_Ycor,
    TEYGr = TEYG-YextG, 
    EYG = 100*EYG_l/YG_l_Ycor,
    EUYG = 100*EUYG_l/YG_l_Ycor,
    YseedG = 100*YseedG_l/YG_l_Ycor,
    YimprG = 100*YimprG_l/YG_l_Ycor,
    EUYGr = EUYG - (YseedG + YimprG),
    TYG = 100*TYG_l/YG_l_Ycor,
    YG = 100*(TEYG_l + EYG_l + EUYG_l + TYG_l)/YG_l_Ycor) %>%
  dplyr::select(ZONE, TEYG:YG) 

# Load SPAM
SPAM <- readRDS(file.path(root, "Cache/SPAMData_ETH.rds"))  %>%
  select(ZONE = zone, Y_SPAM = yield, PROD = TargetProduction) %>%
  mutate(ZONE = toupper(ZONE),
         ZONE = dplyr::recode(ZONE, "BENISHANGUL GUMUZ" = "BENSHANGULGUMUZ",
                                    "SOMALIE" = "SOMALI")) %>%
  filter(ZONE != "ADDIS ABEBA")

# Calculate potentia production per zone
tot_add_prod <- yld_level %>%
  filter(ZONE != "Total") %>%
  select(ZONE, PY) %>%
  left_join(.,SPAM) %>%
  mutate(tot_add_prod = (PY/Y_SPAM*PROD)-PROD) %>%
  select(ZONE, tot_add_prod)

# Calculate production
prod <- SPAM %>%
  select(ZONE, add_prod = PROD) %>%
  mutate(yg_share = "PROD")

# Closing of yield gap assuming decomposition of levels.
yg_close <- yg_share %>%
  filter(ZONE != "Total") %>%
  gather(yg_share, yg, -ZONE) %>%
  left_join(.,tot_add_prod) %>%
  mutate(add_prod = tot_add_prod*yg/100) %>%
  select(ZONE, yg_share, add_prod) %>%
  bind_rows(.,prod) %>%
  group_by(yg_share) %>%
  summarize(yg_close = sum(add_prod)/1000000)

write_csv(yg_close, "Cache/yg_close.csv")

GapClose3a <- GapClose3 %>% 
  summarize(PROD = sum(PROD/1000000, na.rm = TRUE), # in million tons
            TEYG_close = sum(TEYG_close/1000000, na.rm = TRUE),
            EYG_close = sum(EYG_close/1000000, na.rm = TRUE),
            TYG_close = sum(TYG_close/1000000, na.rm = TRUE),
            EUYG_close = sum(EUYG_close/1000000, na.rm = TRUE), 
            POTPROD = sum(POTPROD/1000000, na.rm = TRUE)) 



# Merge SPAM and yield levels
yld_levels <- db4 %>% 
    select(ZONE, Y, Ycor, Yext, TEY, EY, Yseed, Yimpr, PFY, PY, area) %>%
    gather(yldlevel, yld, -ZONE, -area) %>%
    group_by(yldlevel, ZONE) %>%
    summarize(yld = (sum((yld)*area)/sum(area)))

# potential yield
py <- db4 %>% 
  select(ZONE, PY, Ycor, area) %>%
  group_by(ZONE) %>%
  summarize(PY = (sum((PY)*area)/sum(area)),
            Ycor = (sum((Ycor)*area)/sum(area)))

  
rebase <- left_join(py, SPAM) %>%
  ungroup() %>%
  mutate(factor = (PY-Ycor)/(PY-Y_SPAM)) %>%
  select(ZONE, PY, Ycor, Y_SPAM, factor)

rebase2 <- left_join(rebase, yld_levels) %>%
  mutate(yld2 = ((yld-Ycor)/factor)+Y_SPAM) %>%
  select(ZONE, yldlevel, yld2) %>%
  spread(yldlevel, yld2) %>%
  mutate(
    TEYG_l = TEY-Ycor,     # Technical efficiency yield gap using Ycor as basis
    TEYG_s = Ycor/TEY,     # Technical efficiency yield gap using Ycor as basis
    YextG_l = Yext-Ycor,
    EYG_l = EY-TEY,        # Economic yield gap
    EYG_s = TEY/EY,        # Economic yield gap
    EUYG_l = PFY-EY,       # Feasible yield gap
    EUYG_s = EY/PFY,       # Feasible yield gap
    YseedG_l = Yseed-EY,
    YimprG_l = Yimpr-Yseed,
    TYG_l = PY-PFY,        # Technology yield gap
    TYG_s = PFY/PY,        # Technology yield gap
    YG_l = PY-Y,           # Yield gap
    YG_s = Y/PY,           # Yield gap
    YG_l_Ycor = PY-Ycor,   # Yield gap with Ycor as reference
    YG_s_Ycor = Ycor/PY)   # Yield gap with Ycor as reference
    



# Table with absolute yield gap information per zone
# Note that by definition, YG_s computed by weighting individual YG_s values is not the same as multiplication of weighted TEYG_s etc.
# We therefore calculate YG_s as the product of the weighted components.
yg_l <- db4 %>% 
    select(ZONE, TEYG_l, YextG_l, EYG_l, YseedG_l, EUYG_l, TYG_l, YG_l_Ycor, area) %>%
    gather(yg_type, yg, -ZONE, -area) %>%
    group_by(yg_type, ZONE) %>%
    summarize(yg = (sum((yg)*area)/sum(area))) %>%
    ungroup() %>%
    group_by %>%
    mutate(yg_sh = yg/yg[yg_type == "YG_l_Ycor"]) %>%
    select(-yg)


check <- left_join(py, SPAMData) %>%
  mutate(ADDPROD = (PY/Y_SPAM * PROD)-PROD)
        

check2 <- left_join(yg_l, check) %>%
  ungroup() %>%
  mutate(yg_close = ADDPROD*yg_sh) %>%
  group_by(yg_type) %>%
  summarize(yg_close = sum(yg_close)/1000000)
  

%>%
  group_by(yld_level) %>%
  summarize(POTPROD= sum(POTPROD)/1000000)

%>%
    spread(yg_type, yg_sh)

check <- left_join(yg_l, SPAMData) %>%
  
  mutate(POTPROD = rel_yld * PROD) %>%
  group_by(yld_level) %>%
  summarize(POTPROD= sum(POTPROD)/1000000)

  
  
  group_by(Zone) %>%
    summarize(ERROR_l =(sum((ERROR_l)*area, na.rm=TRUE)/sum(area)),
              TEYG_l = (sum((TEYG_l)*area, na.rm=TRUE)/sum(area)),
              EYG_l = (sum((EYG_l)*area, na.rm=TRUE)/sum(area)),
              EUYG_l = (sum((EUYG_l)*area, na.rm=TRUE)/sum(area)),
              TYG_l = (sum((TYG_l)*area, na.rm=TRUE)/sum(area)),
              YG_l = (sum((YG_l)*area, na.rm=TRUE)/sum(area)),
              YG_l_Ycor = (sum((YG_l_Ycor)*area, na.rm=TRUE)/sum(area)),
              YG_lcheck = (ERROR_l+TEYG_l+EYG_l+EUYG_l+TYG_l)),
  db4 %>% 
    dplyr::select(Zone = ZONE, ERROR_l, TEYG_l, EYG_l, EUYG_l, TYG_l, YG_l_Ycor, YG_l, area) %>%
    summarize(Zone = "Total", 
              ERROR_l =(sum((ERROR_l)*area, na.rm=TRUE)/sum(area)),
              TEYG_l = (sum((TEYG_l)*area, na.rm=TRUE)/sum(area)),
              EYG_l = (sum((EYG_l)*area, na.rm=TRUE)/sum(area)),
              EUYG_l = (sum((EUYG_l)*area, na.rm=TRUE)/sum(area)),
              TYG_l = (sum((TYG_l)*area, na.rm=TRUE)/sum(area)),
              YG_l = (sum((YG_l)*area, na.rm=TRUE)/sum(area)),
              YG_l_Ycor = (sum((YG_l_Ycor)*area, na.rm=TRUE)/sum(area)),
              YG_lcheck = (ERROR_l+TEYG_l+EYG_l+EUYG_l+TYG_l))) %>%
  dplyr::select(-ERROR_l,-YG_l, -YG_lcheck)


ZonalYieldGap_l_sh <- ZonalYieldGap_l %>%
  mutate(
    TEYG = 100*TEYG_l/YG_l_Ycor,
    EYG = 100*EYG_l/YG_l_Ycor,
    EUYG = 100*EUYG_l/YG_l_Ycor,
    TYG = 100*TYG_l/YG_l_Ycor,
    YG = 100*(TEYG_l + EYG_l + EUYG_l + TYG_l)/YG_l_Ycor) %>%
  dplyr::select(-TEYG_l:-YG_l_Ycor) 




# http://www.r-bloggers.com/waterfall-plots-in-r/
# Create database for waterfall plot
# Add error, which is very small to target production
wf.df <- GapClose3a %>% 
  dplyr::select(PROD,  TEYG_close, EYG_close, EUYG_close, TYG_close, POTPROD)

wf.df <- as.data.frame(t(wf.df)) %>%
  mutate(category =c("Actual \n production", " Closing \n technical efficiency \n yield gap", " Closing \n allocative \n yield gap",
                     " Closing \n economic \n yield gap", " Closing \n technical \n yield gap", "Potential \n production"),
         sector = category) %>%
  rename(value = V1)

# Create waterfall plot
cbPalette <- c("#009E73", "#CC79A7", "#0072B2", "#D55E00", "black", "#999999")
#waterfall_f(wf.df)

## Determines the spacing between columns in the waterfall chart
offset <- 0.3

Fig_waterfall <- waterfall_f(wf.df, offset=offset) +
  scale_fill_manual(guide="none", values=cbPalette)+
  labs(x="", y="Maize production (million tons)") +
  scale_y_continuous(breaks=seq(0, 45, 5), labels = comma) +
  theme_classic() 




