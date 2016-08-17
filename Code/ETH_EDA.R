###################################################
### EXPLORATORY ANALYSIS of ETH PANEL DATA data ###
###################################################

# From the BID: Addis Ababa, Amhara, Oromiya, SNNP, Tigray, and "other regions" (including Dire Dawa). 
# In those regions with the greatest urban populations, Addis Ababa and Oromiya, 20 EAs were selected; while in all other strata, 15 EAs were selected.

library(plyr)
library(dplyr)
library(ggplot2)
library(stargazer)
library(haven)
library(tidyr)
library(xtable)
library(DescTools)
library(sjPlot)
library(xda)
library(car)

options(scipen=999)

#######################################
####### EXPLORATORY ANALYSIS###########
#######################################

db0 <- dbP
# summary statistics
stargazer(as.data.frame(db0), type = "text") # as.data.frame needed because file are tbl format.
numSummary(db0)
charSummary(db0)

ETH2011 <- filter(db0, surveyyear == 2011)
ETH2013 <- filter(db0, surveyyear == 2013)
summary(ETH2011)
summary(ETH2013)

# Number of observation per hhid in sample
plotobs <- db0 %>%
  group_by(hhid) %>%
  summarize(n=n())
table(plotobs$n)
sum(plotobs$n)

# Check individual variables.
# Number of plots per hh per year
plotnumber <- db0 %>%
            group_by(hhid) %>%
            summarize(n=n())

# Fung and herb
summary(ETH2011$fung)
summary(ETH2013$herb)
table(ETH2013$fung) # Many missing values, not usefull
table(ETH2013$herb) # Many missing values, not usefull

summary(ETH2013$fung)
summary(ETH2011$herb)
table(ETH2011$fung) # Many missing values, not usefull
table(ETH2011$herb) # Many missing values, not usefull

# Yield
sjp.frq(db0$logyld, 
        type = "dens", 
        normal.curve = TRUE, 
        normal.curve.alpha = .3)
Freq(db0$yld)

# Area
sjp.frq(db0$logarea, 
        type = "dens", 
        normal.curve = TRUE, 
        normal.curve.alpha = .3)
Freq(db0$area_gps)

# Labour
sjp.frq(db0$loglab, 
        type = "dens", 
        normal.curve = TRUE, 
        normal.curve.alpha = .3)

sjp.frq(db0$logae, 
        type = "dens", 
        normal.curve = TRUE, 
        normal.curve.alpha = .3)

sjp.frq(db0$hirelab_sh, 
        type = "dens", 
        normal.curve = TRUE, 
        normal.curve.alpha = .3)
Freq(db0$hirelab_sh)

# Fertilizer
sjp.frq(db0$logN, 
        type = "dens", 
        normal.curve = TRUE, 
        normal.curve.alpha = .3)
Freq(db0$N)

# Seed
sjp.frq(db0$seed_q, 
        type = "dens", 
        normal.curve = TRUE, 
        normal.curve.alpha = .3)
Freq(db0$seed_q)
sjp.frq(db0$logseed_q, 
        type = "dens", 
        normal.curve = TRUE, 
        normal.curve.alpha = .3)

# Compare extension and extensions2
xtabs(~db0$extension + db0$extension2)


# Infrastructure
Freq(db0$infra_dummy_finance_ph)
infra <- db0 %>% group_by(ZONE) %>% summarize(infra = mean(infra_dummy_finance_ph),
                                              n=n())

# Visual relationship between variables
Plot(dbP[c("yld", "N")], "yld")
bivariate(dbP, "yld", "N")
names(dbP)

# Scatter point and smoothed line
ggplot(data = filter(db0, N>0), aes(x = N, y = yld)) + geom_point() + geom_smooth() + facet_wrap(~surveyyear)
ggplot(data = db0, aes(x = logN, y = logyld)) + geom_point() + geom_smooth()
ggplot(data = filter(db0, surveyyear == 2012), aes(x = seed_q, y = yld)) + geom_point() + geom_smooth()
ggplot(data = filter(db0), aes(x = N, y = rain_wq)) + geom_point() + geom_smooth() + facet_wrap(~surveyyear, scales = "free")


# Boxplots
ggplot(data = db0, aes(x = factor(seed), y = yld)) + geom_boxplot()
ggplot(data = db0, aes(x = factor(mech), y = yld)) + geom_boxplot()
ggplot(data = db0, aes(x = factor(antrac), y = yld)) + geom_boxplot()
ggplot(data = dbP, aes(x = factor(ext_dummy_pp), y = yld)) + geom_boxplot() + facet_wrap(~surveyyear)
ggplot(data = dbP, aes(x = factor(ext_dummy_ph), y = yld)) + geom_boxplot() + facet_wrap(~surveyyear)
ggplot(data = db0, aes(x = factor(com_type_agricoop_ph), y = yld)) + geom_boxplot() + facet_wrap(~surveyyear)
ggplot(data = db0, aes(x = factor(infra_dummy_finance_ph), y = yld)) + geom_boxplot() + facet_wrap(~surveyyear)

# Cross tables for dichotomous variables
PercTable(db0$mech, db0$antrac,  margins=c(1,2),  rfrq="110", freq=T) 
PercTable(db0$bank_account_own_pp, db0$infra_dummy_finance_ph,  margins=c(1,2),  rfrq="110", freq=T) 
PercTable(db0$ZONE, db0$inter_crop,  margins=c(1,2),  rfrq="110", freq=T) 
PercTable(db0$herb, db0$pest,  margins=c(1,2),  rfrq="110", freq=T) 


