#################################################
####### FIGURES AND TABLES FOR ETHYG ############
#################################################


# PACKAGES ETC 

wdPath <- "D:/Data/Projects/ETHYG/"
setwd(wdPath)

library(plyr)
library(dplyr)
library(ggplot2)
library(stargazer)
library(haven)
library(tidyr)
library(xtable)

options(scipen=999)

# LOAD DATA
# Can be sourced later
db1 <- readRDS("Cache/db1.rds")
db9 <- readRDS("Cache/db9.rds")
db_sfaCD_CRE_Z <- readRDS("Cache/db_sfaCD_CRE_Z.rds")
Prices <- readRDS("Cache/Prices_ETH.rds") %>% rename(Pm = maize, Pn = fertilizer)
source("Code/waterfall_plot.R")
source("Code/sfaTable.R")

# SUMMARY STATISTICS

# Number of unique households
length(unique(db9$hhid))

# summary statistics
dbsum <- db_sfaCD_CRE_Z %>% dplyr::select(Yield = yld, ImprovedSeeds = impr, Slope = slope, 
                                      yesN, Nitrogen = N, Rain = rain_wq, Irrigation = irrig, SOC = SOC2, phdum, 
                                      Labour = lab,  Area = area, crop_count2, AEZ, surveyyear2) 

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

# NOTE THIS IS THE FULL DATABASE NOT WHAT WE USED
Pricessum <- bind_rows(
  Prices %>% 
    dplyr::select(ZONE, Pn, Pm) %>% 
    do(filter(., complete.cases(.))) %>%
    group_by(ZONE) %>%
    summarize(
      NitrogenPrice = mean(Pn, na.rm=T),
      MaizePrice = round(mean(Pm, na.rm=T), digits=0)),
  Prices %>% 
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
Zonalsum <- xtable(Zonalsum, digits = c(0,0,0,0,0,0,0,0,0))
print(Zonalsum, type="html", file=".\\FigTab\\Zonal.html")


# SFA table
sfaCD_CRE_Z <- sfa(logyld ~ noN + logN + loglab + dumoxen + 
                     logarea + irrig + impr + slope + elevation + SOC2 + phdum2 + 
                     rain_wq + AEZ +
                     crop_count2 + surveyyear2 + 
                     noN_bar + logN_bar + loglab_bar + logarea_bar + oxen_bar +
                     irrig_bar + impr_bar + slope_bar + rain_wq_bar + crop_count_bar + 
                     r | 
                     sex + age +  title + literate + ed_any + extension + credit +
                     dist_hh + dist_market + popEA -1
                   ,data = db1, maxit = 1500, restartMax = 20, tol = 0.000001)
summary(sfaCD_CRE_Z, extraPar = TRUE)
lrtest(sfaCD_CRE_Z)

xtable <- sfaTable_f(sfaCD_CRE_Z)[[1]]
ztable <- sfaTable_f(sfaCD_CRE_Z)[[2]]


# 
# Table with yield levels
YieldLevels <- bind_rows(
  db9 %>% 
    dplyr::select(Zone = ZONE, Y, Ycor, TEY, EY, PFY, PY, area) %>%
    group_by(Zone) %>%
    summarize(Y =(sum((Y)*area)/sum(area)),
              Ycor = (sum((Ycor)*area)/sum(area)),
              TEY = (sum((TEY)*area)/sum(area)),
              EY = (sum((EY)*area)/sum(area)),
              PFY = (sum((PFY)*area)/sum(area)),
              PY = (sum((PY)*area)/sum(area))
    ),
  db9 %>% 
    dplyr::select(Zone = ZONE, Y, Ycor, TEY, EY, PFY, PY, area) %>%
    summarize(Zone = "Total", 
              Y =(sum((Y)*area)/sum(area)),
              Ycor = (sum((Ycor)*area)/sum(area)),
              TEY = (sum((TEY)*area)/sum(area)),
              EY = (sum((EY)*area)/sum(area)),
              PFY = (sum((PFY)*area)/sum(area)),
              PY = (sum((PY)*area)/sum(area)))) %>%
  dplyr::select(Zone, Y, Ycor, TEY, EY, PFY, PY)


YieldLevels <- xtable(YieldLevels, digits = c(0,0,0,0,0,0,0,0))
print(YieldLevels, type="html", file=".\\FigTab\\YieldLevels.html")


# Table with relative yield gap information per zone
# Note that by definition, YG_s computed by weighting individual YG_s values is not the same as multiplication of weighted TEYG_s etc.
# We therefore calculate YG_s as the product of the weighted components.
ZonalYieldGap_s <- bind_rows(
  db9 %>% 
    dplyr::select(Zone = ZONE, ERROR_s, TEYG_s, EYG_s, EUYG_s, TYG_s, YG_s_Ycor, YG_s, area) %>%
    group_by(Zone) %>%
    summarize(ERROR_s =(sum((ERROR_s)*area)/sum(area)),
              TEYG_s = (sum((TEYG_s)*area)/sum(area)),
              EYG_s = (sum((EYG_s)*area)/sum(area)),
              EUYG_s = (sum((EUYG_s)*area)/sum(area)),
              TYG_s = (sum((TYG_s)*area)/sum(area)),
              YG_s = (sum((YG_s)*area)/sum(area)),
              YG_s_Ycor = (sum((YG_s_Ycor)*area)/sum(area)),
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
              ERROR_s =(sum((ERROR_s)*area)/sum(area)),
              TEYG_s = (sum((TEYG_s)*area)/sum(area)),
              EYG_s = (sum((EYG_s)*area)/sum(area)),
              EUYG_s = (sum((EUYG_s)*area)/sum(area)),
              TYG_s = (sum((TYG_s)*area)/sum(area)),
              YG_s = (sum((YG_s)*area)/sum(area)),
              YG_s_Ycor = (sum((YG_s_Ycor)*area)/sum(area)),
              ERROR = (1-ERROR_s)*100,
              TEYG = (1-TEYG_s)*100,
              EYG = (1-EYG_s)*100,
              EUYG = (1-EUYG_s)*100,
              TYG = (1-TYG_s)*100,
              YG = (1-(ERROR_s*TEYG_s*EYG_s*EUYG_s*TYG_s))*100,
              YG_Ycor = (1-(TEYG_s*EYG_s*EUYG_s*TYG_s))*100)) %>%
  dplyr::select(Zone, TEYG, EYG, EUYG, TYG, YG = YG_Ycor)


ZonalYieldGap_s <- xtable(ZonalYieldGap_s, digits = c(0,0,0,0,0,0,0))
print(ZonalYieldGap_s, type="html", file=".\\FigTab\\ZonalYG_s.html")

# Table with absolute yield gap information per zone
# Note that by definition, YG_s computed by weighting individual YG_s values is not the same as multiplication of weighted TEYG_s etc.
# We therefore calculate YG_s as the product of the weighted components.
ZonalYieldGap_l <- bind_rows(
  db9 %>% 
    dplyr::select(Zone = ZONE, ERROR_l, TEYG_l, EYG_l, EUYG_l, TYG_l, YG_l_Ycor, YG_l, area) %>%
    group_by(Zone) %>%
    summarize(ERROR_l =(sum((ERROR_l)*area)/sum(area)),
              TEYG_l = (sum((TEYG_l)*area)/sum(area)),
              EYG_l = (sum((EYG_l)*area)/sum(area)),
              EUYG_l = (sum((EUYG_l)*area)/sum(area)),
              TYG_l = (sum((TYG_l)*area)/sum(area)),
              YG_l = (sum((YG_l)*area)/sum(area)),
              YG_l_Ycor = (sum((YG_l_Ycor)*area)/sum(area)),
              YG_lcheck = (ERROR_l+TEYG_l+EYG_l+EUYG_l+TYG_l)),
  db9 %>% 
    dplyr::select(Zone = ZONE, ERROR_l, TEYG_l, EYG_l, EUYG_l, TYG_l, YG_l_Ycor, YG_l, area) %>%
    summarize(Zone = "Total", 
              ERROR_l =(sum((ERROR_l)*area)/sum(area)),
              TEYG_l = (sum((TEYG_l)*area)/sum(area)),
              EYG_l = (sum((EYG_l)*area)/sum(area)),
              EUYG_l = (sum((EUYG_l)*area)/sum(area)),
              TYG_l = (sum((TYG_l)*area)/sum(area)),
              YG_l = (sum((YG_l)*area)/sum(area)),
              YG_l_Ycor = (sum((YG_l_Ycor)*area)/sum(area)),
              YG_lcheck = (ERROR_l+TEYG_l+EYG_l+EUYG_l+TYG_l))) %>%
  dplyr::select(-ERROR_l,-YG_l, -YG_lcheck)

ZonalYieldGap_l <- xtable(ZonalYieldGap_l, digits = c(0,0,0,0,0,0,0))
print(ZonalYieldGap_l, type="html", file=".\\FigTab\\ZonalYG_l.html")

ZonalYieldGap_l_sh <- ZonalYieldGap_l %>%
  mutate(
    TEYG = 100*TEYG_l/YG_l_Ycor,
    EYG = 100*EYG_l/YG_l_Ycor,
    EUYG = 100*EUYG_l/YG_l_Ycor,
    TYG = 100*TYG_l/YG_l_Ycor,
    YG = 100*(TEYG_l + EYG_l + EUYG_l + TYG_l)/YG_l_Ycor) %>%
  dplyr::select(-TEYG_l:-YG_l_Ycor)


ZonalYieldGap_l_sh <- xtable(ZonalYieldGap_l_sh, digits = c(0,0,0,0,0,0,0))
print(ZonalYieldGap_l_sh, type="html", file=".\\FigTab\\ZonalYG_l_sh.html")

# Calculation of potential increase in production when gap is closed on the basis of sample
GapClose1 <- mutate(db9, PROD = Y * area,
                    ERROR_close = ERROR_l*area,
                    TEYG_close = TEYG_l*area,
                    EYG_close = EYG_l*area,
                    EUYG_close = EUYG_l*area,
                    TYG_close = TYG_l*area,
                    POTPROD = PROD + ERROR_close + TEYG_close + EYG_close + TYG_close + EUYG_close)

# Total increase in yield per year
GapClose1a <- GapClose1 %>%
  #group_by(ZONE) %>%
  summarize(PROD = sum(PROD, na.rm=T)/1000000,
            ERROR_close = sum(ERROR_close, na.rm=T)/1000000,
            TEYG_close = sum(TEYG_close, na.rm=T)/1000000,
            EYG_close = sum(EYG_close, na.rm=T)/1000000,
            EUYG_close = sum(EUYG_close, na.rm=T)/1000000,
            TYG_close = sum(TYG_close, na.rm=T)/1000000,
            POTPROD = sum(POTPROD, na.rm=T)/1000000)

# Calculation of potential increase in production when gap is closed on the basis of SPAM and FAO data weighted over two surveyyears.
# As the yield in the LSMS is much lower than the FAO/SPAM yield and also that of the LSMS we apply the different yield shares as found above to the base yield of SPAM.
# We use Ycor as base. This means we assume there was an error (e) and this is corrected for.

SPAMData <- readRDS("./Cache/SPAMData_ETH.rds")  %>%
  rename(ZONE = zone, Y_SPAM = yield, PROD = TargetProduction) %>%
  mutate(ZONE = toupper(ZONE))


# Closing of yield gaps per zone
# Note that for some regions, notably those with very low potential yield (Central and Western), closing TEY and EY already results in
# Closing the gap. To avoid negative closing the EYG, EUYG and TYG are capped.
# The reason for overshooting can be caused by a variety of factors, including mismeasurement. Most likely is that Nyp is too high for regions
# With a very low potential. The all over impact is low as the involved regions have very limited maize production.


GapClose2 <- db9 %>% 
  group_by(ZONE) %>%
  summarize(
    TEYG_s = sum(TEYG_s*area)/sum(area),
    EYG_s = sum(EYG_s*area)/sum(area),
    #TYG_s = (sum((TYG_s)*area)/sum(area)), # TYG_s based on LSMS yield, not used
    #YG_s = (sum((YG_s)*area)/sum(area)), # YG_s based on LSMS yield, not used
    EUYG_s = sum(EUYG_s*area)/sum(area),
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
  summarize(PROD = sum(PROD/1000000), # in million tons
            TEYG_close = sum(TEYG_close/1000000),
            EYG_close = sum(EYG_close/1000000),
            TYG_close = sum(TYG_close/1000000),
            EUYG_close = sum(EUYG_close/1000000), 
            POTPROD = sum(POTPROD/1000000)) %>%
  mutate(check2 = TEYG_close + EYG_close + EUYG_close + TYG_close+PROD)

# Closing of yield gap assuming decomposition of levels.

GapClose3 <- SPAMData %>%
  rename(Zone = ZONE) %>%
  left_join(ZonalYieldGap_l_sh, .) %>%
  left_join(.,YieldLevels) %>%
  filter(Zone != "Total") %>%
  mutate(POTPROD = PY/Y_SPAM*PROD, # Total production when YG is closed
         TEYG_close = (POTPROD-PROD)*TEYG/100, # Additional production when closing TEYG
         EYG_close = (POTPROD-PROD)*EYG/100, # Additional production when closing EYG
         EUYG_close = (POTPROD-PROD)*EUYG/100, # Additional production when closing EUYG
         TYG_close = (POTPROD-PROD)*TYG/100, # Additional production when closing TYG
         check2 = TEYG_close + EYG_close + EUYG_close + TYG_close+PROD)

GapClose3a <- GapClose3 %>% 
  summarize(PROD = sum(PROD/1000000, na.rm = TRUE), # in million tons
            TEYG_close = sum(TEYG_close/1000000, na.rm = TRUE),
            EYG_close = sum(EYG_close/1000000, na.rm = TRUE),
            TYG_close = sum(TYG_close/1000000, na.rm = TRUE),
            EUYG_close = sum(EUYG_close/1000000, na.rm = TRUE), 
            POTPROD = sum(POTPROD/1000000, na.rm = TRUE)) 

# http://www.r-bloggers.com/waterfall-plots-in-r/
# Create database for waterfall plot
# Add error, which is very small to target production
wf.df <- GapClose3a %>% 
  dplyr::select(PROD,  TEYG_close, EYG_close, EUYG_close, TYG_close, POTPROD)

wf.df <- as.data.frame(t(wf.df)) %>%
  mutate(category =c("Actual \n production", " Closing \n technical efficiency \n yield gap", " Closing \n economic \n yield gap",
                     " Closing \n feasible \n yield gap", " Closing \n technical \n yield gap", "Potential \n production"),
         sector = category) %>%
  rename(value = V1)

# Create waterfall plot
cbPalette <- c("#009E73", "#CC79A7", "#0072B2", "#D55E00", "black", "#999999")
#waterfall_f(wf.df)

## Determines the spacing between columns in the waterfall chart
offset <- 0.3

waterfall <- waterfall_f(wf.df, offset=offset) +
  scale_fill_manual(guide="none", values=cbPalette)+
  labs(x="", y="Maize production (million tons)") +
  scale_y_continuous(breaks=seq(0, 40, 5), labels = comma) +
  theme_classic() 


#print(waterfall)
library(Cairo)
ggsave(plot = waterfall, ".\\FigTab\\ETH_Waterfall.png", height = 150, width = 200, type = "cairo-png", units="mm")


# Distribution of relative yield gaps
db11 <- db9%>%
  dplyr::select(ZONE, REGNAME, surveyyear, ERROR_s, TEYG_s, EYG_s, EUYG_s, TYG_s, YG_s_Ycor) %>%
  gather(yieldgap, value, ERROR_s:YG_s_Ycor) %>%
  filter(yieldgap!="ERROR_s") %>% 
  droplevels() %>%
  mutate(value = (1-value)*100,
         yieldgap = factor(yieldgap, levels = c("TEYG_s", "EYG_s", "EUYG_s", "TYG_s",  "YG_s_Ycor")))

# 
# # Note that average error=0 and therefore not interesting to show.
# # Show plot with minimum hinges and no outliers.
# boxplot <-ggplot(data=db11, aes(x=yieldgap, y=value)) +
#   #geom_boxplot(fill=c("#D55E00", "#009E73", "#0072B2", "#CC79A7", "#999999"), outlier.colour = NA) +
#   geom_boxplot(fill=c("#D55E00", "#009E73", "#0072B2", "#CC79A7", "#999999"), outlier.colour = "black") +
#   stat_boxplot(geom ='errorbar') +
#   guides(fill=FALSE) +
#   stat_summary(fun.y=mean, geom="point", shape=5, size=4) +
#   theme_classic() +
#   #geom_jitter(position=position_jitter(width=.1, height=0)) +
#   labs(x="", y="Yield gap (%)")+ 
#   scale_x_discrete(breaks=c("TEYG_s", "EYG_s", "EUYG_s", "TYG_s", "YG_s_Ycor"), 
#                    labels=c("Techical efficiency \n yield gap", "Economic \n yield gap", 
#                             "Feasible \n yield gap", "Technical \n yield gap",  "Total \n yield gap")) 
# #facet_wrap(~zone) +
# #scale_y_continuous(breaks=seq(0, 100, 25)) 
# 
# # Rescale axes.
# #sts <- boxplot.stats(db11$value)$stats  # Compute lower and upper whisker limits
# #boxplot = boxplot + coord_cartesian(ylim = c(-5,max(sts)*1.05))
# boxplot
# ggsave(plot = boxplot, ".\\FigTab\\ETH_Distribution.png", height = 150, width = 200, type = "cairo-png", units="mm")
# 
# # Distribution of absolute yield gaps
# db12 <- db9%>%
#   dplyr::select(ZONE, REGCODE, surveyyear, ERROR_l, TEYG_l, EYG_l, EUYG_l, TYG_l, YG_l_Ycor) %>%
#   gather(yieldgap, value, ERROR_l:YG_l_Ycor) %>%
#   filter(yieldgap!="ERROR_l") %>% 
#   droplevels() %>%
#   mutate(yieldgap = factor(yieldgap, levels = c("TEYG_l", "EYG_l", "TYG_l", "EUYG_l", "YG_l_Ycor")))
# 
# # Note that average error=0 and therefore not interesting to show.
# # Show plot with minimum hinges and no outliers.
# boxplot2 <-ggplot(data=db12, aes(x=yieldgap, y=value)) +
#   geom_boxplot(fill=c("#D55E00", "#009E73", "#0072B2", "#CC79A7", "#999999"), outlier.colour = NA) +
#   #geom_boxplot(fill=c("#D55E00", "#009E73", "#0072B2", "#CC79A7", "#999999"), outlier.colour = "black") +
#   stat_boxplot(geom ='errorbar') +
#   guides(fill=FALSE) +
#   stat_summary(fun.y=mean, geom="point", shape=5, size=4) +
#   theme_classic() +
#   labs(x="", y="Yield gap (%)")+ 
#   scale_x_discrete(breaks=c("TEYG_l", "EYG_l", "TYG_l", "EUYG_l", "YG_l_Ycor"), 
#                    labels=c("Techical efficiency \n yield gap", "Economic \n yield gap", 
#                             "Technical \n yield gap", " Economically \n unexploitable \n yield gap", "Yield gap"))
# #+ scale_y_continuous(breaks=seq(0, 100, 25)) 
# 
# # Rescale axes.
# boxplot2
# ggsave(plot = boxplot2, ".\\FigTab\\ETH_Distribution2.png", height = 150, width = 200, type = "cairo-png", units="mm")
# 

# Maps with GYGA yield potential and plot information
# transform shapefile in dataframe for ggplot. rownames are used as ID for the countries. Using ID gives strange results. 
# Additional data is linked back again
library(rgdal)
library(sp)
library(raster)
library(foreign)
library(gdata)
library(RColorBrewer)

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
meanYield <- ddply(db9,.(lon, lat), summarize, meanYield = (sum(Y*area)/sum(area))/1000)
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

#GYGA_LSMS
ggsave(plot = GYGA_LSMS, ".\\FigTab\\GYGA_LSMS.png", height = 150, width = 200, type = "cairo-png", units="mm")

# Zonal map with community yield levels
# read in map of tanzania as a SpatialPolygonsDataFrame

countryMap <- getData('GADM', country = "ETH", level = 1) 

# Rename zones using LSMS names
countryMap@data$ZONE <- countryMap@data$NAME_1
countryMap@data$ZONE[countryMap@data$NAME_1 %in% c("Oromia")] <- "Oromiya"
countryMap@data$ZONE[countryMap@data$NAME_1 %in% c("Somali")] <- "Somalie"
countryMap@data$ZONE[countryMap@data$NAME_1 %in% c("Benshangul-Gumaz")] <- "Benishangul Gumuz"
countryMap@data$ZONE[countryMap@data$NAME_1 %in% c("Southern Nations, Nationalities and Peoples")] <- "SNNP"
countryMap@data$ZONE[countryMap@data$NAME_1 %in% c("Gambela Peoples")] <- "Gambella"
countryMap@data$ZONE[countryMap@data$NAME_1 %in% c("Harari People")] <- "Harari"
countryMap@data$ZONE <- factor(countryMap@data$ZONE)

# Remove Addis Ababa and Dire Dawa
countryMap <- countryMap[!(countryMap@data$ZONE %in% c("Addis Abeba", "Dire Dawa")),]
plot(countryMap)

#  fortify spatial data to use with ggplot and join using join functions from dplyr
#    The join is on id, make sure all ids are character vectors
tf <- fortify(countryMap)
countryMap@data <- rename(countryMap@data, id = ID_1)
countryMap@data$id <- as.character(countryMap@data$id)
tf2 <- left_join(tf, countryMap@data)

# Use ggplot to plot map of Tanzania, specifying the labels and choosing nice colours
#    from the RColorBrewer package

#display.brewer.all()

ZONE_LSMS <- ggplot()+
  geom_polygon(data=tf2, aes(x=long, y=lat, group=group, fill=ZONE), colour="black")+
  geom_point(data=meanYield, aes(x=lon, y=lat, size=(meanYield2)), colour="black")+
  scale_fill_brewer(name = "Zones", palette = "Set1") +
  scale_size_manual(name="Average yield (tons)", values=c(1.5, 2.5, 3.5, 4,5)) +
  coord_equal()+
  labs(x="", y="")+
  theme_classic() +
  theme(legend.key=element_blank(),
        line = element_blank(),
        axis.text = element_blank())
#ZONE_LSMS 

ggsave(plot = ZONE_LSMS, ".\\Cache\\ZONE_LSMS.png", height = 150, width = 200, type = "cairo-png", units="mm")
save.image("Cache/FigTab.Rdata")


