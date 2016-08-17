# -------------------------------------
# area imputation file for ETHiopia 2013
# missing gps measurements

library(haven)
library(dplyr)

dataPath <- "C:/Users/Tomas/Documents/LEI/data/ETH"

# -------------------------------------
# plot characteristics

plots <- read_dta(file.path(dataPath, "Post-Planting/sect3_pp_w2.dta")) %>%
  dplyr::select(holder_id, household_id2, parcel_id, field_id, status = pp_s3q03,
                irrigated = pp_s3q12, urea = pp_s3q15, dap = pp_s3q18,
                other_inorg = pp_s3q20a, manure = pp_s3q21, compost = pp_s3q23,
                other_org = pp_s3q25)

plots$rented <- ifelse(plots$status %in% 6, 1, 0) # 6 for rentedout, household use otherwise
plots$fallowed <- ifelse(plots$status %in% 3, 1, 0) # 3 if plot left fallow
plots$irrigated <- ifelse(plots$irrigated %in% 1, 1, 0)

# fertilizer - inorganic

plots$dap <- ifelse(plots$dap %in% 1, 1, 0)
plots$urea <- ifelse(plots$urea %in% 1, 1, 0)
plots$other_inorg <- ifelse(plots$other_inorg %in% 1, 1, 0)
plots$inorgFert <- with(plots, ifelse(dap == 1 | urea == 1 | other_inorg ==1, 1, 0))

# fertilizer - organic

plots$manure <- ifelse(plots$manure %in% 1, 1, 0)
plots$compost <- ifelse(plots$compost %in% 1, 1, 0)
plots$other_org <- ifelse(plots$other_org %in% 1, 1, 0)
plots$orgFert <- with(plots,
                      ifelse(manure == 1 | compost == 1 | other_org ==1, 1, 0))

plots <- select(plots, holder_id, household_id2, parcel_id, field_id,
                irrigated, rented, fallowed, inorgFert, orgFert)

# count the number of plots owned by a household

plotCount <- group_by(plots, household_id2) %>% summarise(n=n())

# join with other plot information

plots <- left_join(plots, plotCount)

# soil quality is a parcel level variable

parcel <- read_dta(file.path(dataPath, "Post-Planting/sect2_pp_w2.dta")) %>%
  dplyr::select(holder_id, household_id2, parcel_id,
                owned = pp_s2q04, soil_qlty=pp_s2q15)

parcel$owned <- ifelse(parcel$owned %in% 1, 1, 0)
parcel$soil_qlty <- as_factor(parcel$soil_qlty)
parcel$soil_qlty <- relevel(parcel$soil_qlty, ref = "Poor")

plots <- left_join(plots, parcel)

# Hired labour

pp_lab <- read_dta(file.path(dataPath, "Post-Planting/sect3_pp_w2.dta")) %>%
  dplyr::select(holder_id, household_id2, parcel_id, field_id, pp_s3q28_a:pp_s3q28_h) %>%
  transmute(holder_id, household_id2, parcel_id, field_id,
            hirM=pp_s3q28_a*pp_s3q28_b,
            hirF=pp_s3q28_d*pp_s3q28_e,
            hirC=pp_s3q28_g*pp_s3q28_h
  )

# make all NA values zero
pp_lab[is.na(pp_lab)] <- 0

# sum all labour across a single plot - all measured in days
pp_lab <- transmute(pp_lab, holder_id, household_id2, parcel_id, field_id,
                  hir_lab=hirM + hirF + hirC)

# post harvest hired labour - measured at the crop level

ph_lab <- read_dta(file.path(dataPath, "Post-Harvest/sect10_ph_w2.dta")) %>%
  dplyr::select(holder_id, household_id2, parcel_id, field_id,
                crop=crop_code, ph_s10q01_b:ph_s10q01_h) %>%
  transmute(holder_id, household_id2, parcel_id, field_id, crop,
            hirM=ph_s10q01_b,
            hirF=ph_s10q01_e,
            hirC=ph_s10q01_h
  )

ph_lab[is.na(ph_lab)] <- 0

# sum all labour across a single crop - all measured in days                   
ph_lab <- transmute(ph_lab, holder_id, household_id2, parcel_id, field_id,
                  ph_hir_lab = hirM + hirF + hirC)

# group by crops to find total labour used on
# a field

ph_lab <- group_by(ph_lab, holder_id, household_id2, parcel_id, field_id) %>%
  summarise(hir_lab = sum(ph_hir_lab))

# bind and sumarise the data
hir_lab <- rbind(pp_lab, ph_lab) %>%
  group_by(holder_id, household_id2, parcel_id, field_id) %>%
  summarise(hiredLabour = sum(hir_lab))

# join the labour with the rest of the plot
# information

plots <- left_join(plots, hir_lab)

rm(list=ls()[!ls() %in% c("plots", "dataPath")])                               

# -------------------------------------
# household characteristics

# plot manager characteristics

pmc <- read_dta(file.path(dataPath, "Household/sect1_hh_w2.dta")) %>%
  filter(hh_s1q02 %in% 1) %>% # 1 for head of household
  dplyr::select(household_id2, individual_id, sex=hh_s1q03,
                age=hh_s1q04_a, rural, weight=pw2)

pmc$sex <- ifelse(pmc$sex %in% 2, 1, 0)
pmc$rural <- ifelse(pmc$rural %in% 1, 1, 0)

# education

# education <- read_dta(file.path(dataPath, "household/sect2_hh_w2.dta")) %>%
#   dplyr::select(household_id2, individual_id, educ = hh_s2q05)
# 
# education$educ <- as_factor(education$educ)

# Join the information on the plot manager

# pmc <- left_join(pmc, education)

# count number of people of each age group in
# the household

ag <- read_dta(file.path(dataPath, "Household/sect1_hh_w2.dta")) %>%
  dplyr::select(household_id2, individual_id, age=hh_s1q04_a, sex=hh_s1q03) %>%
  mutate(age_group = cut(age, c(0, 5, 14, 39, 59, 100),
                         include.lowest = TRUE, right = FALSE))
  
ag$sex <- ifelse(ag$sex %in% 2, 1, 0)

# set ages without sex
ag$less5 <- ifelse(ag$age_group == '[0,5)', 1, NA)
ag$g5l15 <- ifelse(ag$age_group == '[5,14)', 1, NA)
ag$above60 <- ifelse(ag$age_group == '[59,100]', 1, NA)

# ages with sex
ag$male15_39 <- ifelse(ag$age_group=='[14,39)' & ag$sex==0, 1, NA)
ag$female15_39 <- ifelse(ag$age_group=='[14,39)' & ag$sex==1, 1, NA)
ag$male40_59 <- ifelse(ag$age_group=='[39,59)' & ag$sex==0, 1, NA)
ag$female40_59 <- ifelse(ag$age_group=='[39,59)' & ag$sex==1, 1, NA)

# summarise to get variables per household.
by_hhid <- group_by(ag, household_id2) %>% summarise(
  less5 = sum(less5, na.rm = TRUE),
  above60 = sum(above60, na.rm = TRUE),
  male15_39 = sum(male15_39, na.rm = TRUE),
  female15_39 = sum(female15_39, na.rm = TRUE),
  male40_59 = sum(male40_59, na.rm = TRUE),
  female40_59 = sum(female40_59, na.rm = TRUE),
  g5l15 = sum(g5l15, na.rm = TRUE))

# join with plot manager information to get the 
# household characteristics (hc)

hc <- left_join(pmc, by_hhid)
rm(pmc, ag, by_hhid)

# ------------------------------------
# areas

areas <- read_dta(file.path(dataPath, "Post-Planting/sect3_pp_w2.dta")) %>%
  dplyr::select(region=saq01, zone=saq02, woreda=saq03,
                holder_id, household_id2, parcel_id,
                field_id, area_sr=pp_s3q02_a,
                local_unit=pp_s3q02_c, area_gps=pp_s3q05_a)

areas$region <- as.numeric(areas$region)
areas$local_unit <- as.numeric(areas$local_unit)

# it is possible to get some idea of 
# a conversion for most of the local units
# the guys we are missing are 7 and 8

table(areas$local_unit)
sum(table(areas$local_unit)[7:8]) # 8282

# unit conversions

unitCon <- read_dta(file.path(dataPath, "ET_local_area_unit_conversion.dta")) %>%
  dplyr::select(region, zone, woreda, local_unit, conversion)

unitCon$region <- as.numeric(unitCon$region)
unitCon$local_unit <- as.numeric(unitCon$local_unit)

# many missing values when joining areas with
# conversion. Only way to reduce this is to use
# average of local conversion. No regional
# differences.

unitCon <- group_by(unitCon, local_unit) %>%
  summarise(conversion = mean(conversion, na.rm=TRUE))

# join areas with converisons

areas <- left_join(areas, unitCon)

# self reported area is our most important
# variable for imputing the missing areas
# so remove all rows where area_sr is missing

areas <- areas[!is.na(areas$area_sr), ]

# if the local unit is square metres or
# hectacres then conversion is easy

areas$conversion <- ifelse(areas$local_unit %in% 1, 1/0.0001,
                           ifelse(areas$local_unit %in% 2, 1, areas$conversion))

# find out how many plots will still be missing
# gps areas after imputation and remove plots
# for which we do not have a conversion

sum(is.na(areas$area_gps) & is.na(areas$conversion)) # 237

areas <- areas[!is.na(areas$conversion), ]

# multiply the self reported are by the 
# conversion. further multiply by 0.0001
# to get square metres to hectacres

sm2ha <- 0.0001
areas$area_sr <- areas$area_sr*areas$conversion*sm2ha                         
areas$area_gps <- areas$area_gps*sm2ha

# make a plot to compare the self reported and the
# gps measurements. Outliers so use 99th quantile
# of areas_gps. These will be winsored anyway

q99gps <- quantile(areas$area_gps, 0.99, na.rm=TRUE)
with(areas[areas$area_sr < q99gps & areas$area_gps < q99gps,],
     plot(area_gps ~ area_sr))
abline(0, 1, col='red', lwd=3)

# winsor self reported areas as they are likely to
# include outliers - see world bank imputation 
# paper

source("C:/Users/Tomas/Documents/LEI/pro-gap/winsor.R")                 
areas$area_sr <- winsor(areas$area_sr, 0.001, lower=FALSE)

# make another comparison plot
with(areas, plot(area_gps ~ area_sr))
areas <- select(areas, -region, -zone, -woreda, -conversion, -local_unit)

rm(list=ls()[!ls() %in% c("plots", "hc", "areas")])

# -------------------------------------
# join all information together

impData <- left_join(areas, plots) %>% left_join(hc)

rm(list=ls()[!ls() %in% "impData"])

# -------------------------------------
# Imputation code

# install.packages('mice')
# install.packages('VIM')

library(mice)
library(VIM)
library(lattice)

# kill of household id and plotnum. These can't go into the 
# imputation

ids <- impData[, c("holder_id","household_id2", "parcel_id","field_id")]
impData$holder_id <- impData$household_id2 <-
  impData$parcel_id <- impData$field_id <-
  impData$individual_id <- NULL

# following functions give an overview of the
# missingness in our data

md.pairs(impData)
md.pattern(impData)
marginplot(impData[, c("area_sr", "area_gps")], col=c('blue', 'red', 'orange'))

# create interaction variables to include in the 
# imputation

impData <- mutate(impData,
                  area_sr2 = area_sr^2,
                  area_sr3 = area_sr^3,
                  area_sr_rural = area_sr*rural,
                  area_sr_sex = area_sr*sex,
                  area_sr_rented = area_sr*rented,
                  area_sr_fallowed = area_sr*fallowed,
                  area_sr_age = area_sr*age,
                  area_sr_owned = area_sr*owned
                  )

# MICE imputation method

imp <- mice(impData, m=50)

# create new areas dataset with new imputed values

imputed <- complete(imp)
areas <- select(impData, area_gps, area_sr)
areas$area_gps_mi50 <- imputed$area_gps

# rejoin with id variables

areas <- cbind(ids, areas)

# save areas to a file

write_dta(areas, "C:/Users/Tomas/Documents/LEI/data/ETH/areas_ETH2013.dta")


