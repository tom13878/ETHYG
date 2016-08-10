#######################################
#### CLEAN DATA PREPARE NGA PANEL #####
#######################################

dataPath <- "N:/Internationaal Beleid  (IB)/Projecten/2285000066 Africa Maize Yield Gap/SurveyData/"
wdPath <- "D:\\Dropbox\\Michiel_research\\2285000066 Africa Maize Yield Gap"
setwd(wdPath)


############################################
############## READ THE DATA ###############
############################################

source("")



#######################################
############## PACKAGES ETC ###########
#######################################

detach(package:dplyr)
library(plyr)
library(dplyr)
library(stargazer)
library(broom)
library(DescTools)
library(ggplot2)
library(xtable)
library(frontier)
library(moments)
library(tidyr)
library(openxlsx)

options(scipen=999)

##################################
############## FIGURES ###########
##################################

