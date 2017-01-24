# -------------------------------------
# Ethiopia data 2011 - 2012 survey
# -------------------------------------

if(Sys.info()["user"] == "Tomas"){
  dataPath <- "C:/Users/Tomas/Documents/LEI/data/ETH/2011/Data"
} else {
  dataPath <- "N:/Internationaal Beleid  (IB)/Projecten/2285000066 Africa Maize Yield Gap/SurveyData/ETH/2011/Data"
}

library(haven)
library(lubridate)
library(dplyr)
library(sjmisc)
library(car)

options(scipen=999)

#######################################
############## LOCATION ###############
#######################################

# file containing the zone, region and district

location <- read_dta(file.path(dataPath, "/sect1_hh_w1.dta")) %>%
  select(household_id, REGCODE = saq01, ZONECODE = saq02, WOREDACODE = saq03, KEBELECODE = saq06, ea_id, rural) 

# recode the rural variable to match with wave 2
location$rural <- ifelse(location$rural %in% 1, 1, 2)
location$type <- factor(location$rural, levels=1:2, labels=c("RURAL", "SMALL TOWN"))
location$rural <- ifelse(location$rural %in% 1, 1, 0)
location$REGCODE <- as.integer(location$REGCODE)
location <- unique(location)

# match up with the names from the survey (prepared in a seperate file)

REGZONE <- read.csv(file.path(paste0(dataPath,"/../../.."), "Other/Spatial/ETH/REGZONEETH.csv"))

# join with household identifications
location <- left_join(location, REGZONE)

rm(REGZONE)

#######################################
########### SOCIO/ECONOMIC ############
#######################################

HH11 <- read_dta(file.path(dataPath, "sect1_hh_w1.dta")) %>%
  select(household_id, individual_id,
         ea_id, status=hh_s1q02, sex=hh_s1q03, age=hh_s1q04_a,
         religion=hh_s1q07, marital=hh_s1q08)

HH11$status <- toupper(as_factor(HH11$status))
HH11$religion <- toupper(as_factor(HH11$religion))
HH11$sex <- recode(HH11$sex, c("2=1; 1=0; else = NA"))
HH11$marital <- toupper(as_factor(HH11$marital)) 

# make a new variable cage (cut age = cage) which splits
# individuals according to their age group with
# breaks at 15, 55 and the max age

HH11$cage <- cut(HH11$age, breaks = c(0, 15, 55, max(HH11$age, na.rm=TRUE)),
                 labels=c("0-15", "16-55", "56+"), include.lowest = TRUE, right = TRUE)

# education of household members and sum
# of education of all household members
# between the ages of 15 and 55

education <- read_dta(file.path(dataPath, "sect2_hh_w1.dta")) %>%
  select(household_id, individual_id, ea_id,
         literate=hh_s2q02, ed_any=hh_s2q03)

education$literate <- ifelse(education$literate %in% 2, 0, education$literate)
education$ed_any <- ifelse(education$ed_any %in% 2, 0, education$ed_any)

# summarise the data: get sum of education
# of household members 15-55 and number of
# household members 15:55

HH11_x <- group_by(HH11, household_id) %>%
  summarise(N1555=sum(cage %in% "16-55"),
            family_size=n())

# -------------------------------------
# death in the family
# -------------------------------------

death <- read_dta(file.path(dataPath, "sect8_hh_w1.dta")) %>%
  select(household_id, code=hh_s8q00, death=hh_s8q01) %>% 
  filter(code %in% c("101", "102")) %>% select(-code) %>%
  group_by(household_id) %>%
  summarise(death=ifelse(any(death %in% 1), 1, 0))

HH11 <- left_join(HH11, education) %>%
  left_join(HH11_x) %>%
  left_join(death); rm(education, HH11_x, death)

#######################################
############### OUTPUT ################
#######################################

# in the 2011 wave of the data, the respondents
# were asked about the output per 2 metre by 2 metre
# area. They were also asked how much the total field
# produced.

oput <- remove_all_labels(read_dta(file.path(dataPath, "sect9_ph_w1.dta"))) %>%
  transmute(household_id, holder_id, parcel_id, field_id,
         crop_code, day_cut=ph_s9q02_a, month_cut=ph_s9q02_b,
         crop_qty_harv_fresh_kg=ph_s9q03_a,
         crop_qty_harv_fresh_g=ph_s9q03_b/1000,
         crop_qty_harv_dry_kg=ph_s9q05_a,
         crop_qty_harv_dry_g=ph_s9q05_b/1000,
         crop_qty_harv_tot_kg=ph_s9q12_a,
         crop_qty_harv_tot_g=ph_s9q12_b/1000,
         harv_month_start=ph_s9q13_a,
         harv_month_end=ph_s9q13_b, crop_name)

oput$crop_code <- as.integer(oput$crop_code)
x <- c("crop_qty_harv_fresh_kg","crop_qty_harv_fresh_g","crop_qty_harv_dry_kg",
       "crop_qty_harv_dry_g","crop_qty_harv_tot_kg","crop_qty_harv_tot_g")
make0 <- is.na(oput[, x])
oput[, x][make0] <- 0

oput <- transmute(oput, household_id, holder_id, parcel_id, field_id,
                  crop_code, day_cut, month_cut, 
                  crop_qty_harv_fresh=crop_qty_harv_fresh_kg+crop_qty_harv_fresh_g,
                  crop_qty_harv_dry=crop_qty_harv_dry_kg+crop_qty_harv_dry_g,
                  crop_qty_harv=crop_qty_harv_tot_kg+crop_qty_harv_tot_g,
                  crop_name)

# -------------------------------------
# add a dummy if a legume was grown
# count number of crops per field
# crop codes taken from appendix in 
# the BID

legumes <- c(11:18, 36, 118)

oput_x <- group_by(oput, holder_id, household_id, parcel_id, field_id) %>%
  summarise(crop_count=sum(!is.na(crop_code)),
            legume = ifelse(any(crop_code %in% legumes), 1, 0))

oput <- left_join(oput, oput_x); rm(oput_x)

# recode observations with quantity 0 to NA again as they have articially set to 0 above
# should also be done for crop cut qty but we do not use them so not yet done.
oput[oput$crop_qty_harv ==0,] <- NA
 
rm("legumes")

# the price of each crop and how much was
# sold is stored in a seperate section of 
# the household (section 11)

oput2 <- read_dta(file.path(dataPath, "sect11_ph_w1.dta")) %>%
  select(holder_id, household_id, crop_code,
         sold=ph_s11q01, sold_qty_kg=ph_s11q03_a, sold_qty_gr=ph_s11q03_b,
         value=ph_s11q04_a, sold_month=ph_s11q06_a, sold_year=ph_s11q06_b,
         trans_cost=ph_s11q09) 

oput2$sold_qty_kg[is.na(oput2$sold_qty_kg)] <- 0
oput2$sold_qty_gr[is.na(oput2$sold_qty_gr)] <- 0
oput2$sold_qty <- oput2$sold_qty_kg + (oput2$sold_qty_gr/1000)
oput2$sold_qty[oput2$sold_qty==0] <- NA
oput2$sold_qty_kg <-NULL
oput2$sold_qty_gr <-NULL
oput2$crop_price <- oput2$value/oput2$sold_qty
oput2$crop_code <- as.integer(oput2$crop_code)
oput2$sold <- toupper(as_factor(oput2$sold))
oput2$sold_month <- month(oput2$sold_month, label=TRUE)

oput <- left_join(oput, oput2) %>% unique; rm(oput2)

#######################################
############## CHEMICAL ###############
#######################################

# -------------------------------------
# variables recorded at parcel, field,
# and crop level

# parcel level
parcel <- remove_all_labels(read_dta(file.path(dataPath, "sect2_pp_w1.dta"))) %>%
  select(household_id, holder_id, parcel_id, fields=pp_s2q02,
                title=pp_s2q04)

parcel$title <- ifelse(parcel$title %in% 2, 0, parcel$title)
parcel$parcel_id <- as.integer(parcel$parcel_id)

# field level variables

field <- read_dta(file.path(dataPath, "sect3_pp_w1.dta")) %>%
  select(holder_id, household_id, parcel_id, field_id,
                crop_stand=pp_s3q10,
                extension=pp_s3q11, irrig=pp_s3q12, fert_any=pp_s3q14,
                manure=pp_s3q21, compost=pp_s3q23,
                other_org=pp_s3q25) 

field$irrig <- ifelse(field$irrig %in% 2, 0, field$irrig)
field$crop_stand <- toupper(as_factor(field$crop_stand))
field$extension <- recode(field$extension, c("1=1; 2=0; 3=NA")) # for some reason there are a few '3' in the variable => NA.
field$parcel_id <- as.integer(field$parcel_id)
field$field_id <- as.integer(field$field_id)

# crop level variables

crop <- read_dta(file.path(dataPath, "sect4_pp_w1.dta")) %>%
  select(household_id, holder_id, parcel_id, field_id, crop_code,
         cropping=pp_s4q02, month=pp_s4q12_a, crop_area=pp_s4q03,
         pest=pp_s4q05, herb=pp_s4q06, fung=pp_s4q07, impr=pp_s4q11,
        crop_name)

crop$pest <- ifelse(crop$pest %in% 2, 0, crop$pest)
crop$herb <- ifelse(crop$herb %in% 2, 0, crop$herb)
crop$fung <- ifelse(crop$fung %in% 2, 0, crop$fung)
crop$impr <- recode(crop$impr, c("2=1; 1=0; 3=NA"))
crop$cropping <- as_factor(crop$cropping)
crop$month <- toupper(as_factor(crop$month))
crop$crop_code <- as.integer(crop$crop_code)
crop <- remove_all_labels(crop)

seed <- remove_all_labels(read_dta(file.path(dataPath, "sect5_pp_w1.dta"))) %>%
  select(household_id, holder_id, parcel_id, field_id, crop_code,
         impr2 = pp_s5q01, seed_q_kilo = pp_s5q19_a, seed_q_gram = pp_s5q19_b) 

seed$seed_q_kilo[is.na(seed$seed_q_kilo)] <- 0
seed$seed_q_gram[is.na(seed$seed_q_gram)] <- 0
seed$seed_q = seed$seed_q_kilo + (seed$seed_q_gram/1000)
seed$seed_q[seed$seed_q == 0] <- NA
seed$seed_q_kilo <-NULL
seed$seed_q_gram <-NULL

# -------------------------------------
# unit of observation is not fertilizer

fert1 <- read_dta(file.path(dataPath, "sect3_pp_w1.dta")) %>%
  transmute(holder_id, household_id, parcel_id, field_id, typ=pp_s3q15,
            qty1=pp_s3q16_a, qty2=pp_s3q16_b/1000)
x <- c("qty1", "qty2")
make0 <- is.na(fert1[, x])
fert1[, x][make0] <- 0
fert1 <- transmute(fert1, holder_id, household_id, parcel_id, field_id, typ,
                   qty=qty1+qty2)
fert1$typ <- ifelse(fert1$typ %in% 1, "UREA", NA)

fert2 <- read_dta(file.path(dataPath, "sect3_pp_w1.dta")) %>%
  transmute(holder_id, household_id, parcel_id, field_id, typ=pp_s3q18,
            qty1=pp_s3q19_a, qty2=pp_s3q19_b/1000)
x <- c("qty1", "qty2")
make0 <- is.na(fert2[, x])
fert2[, x][make0] <- 0
fert2 <- transmute(fert2, holder_id, household_id, parcel_id, field_id, typ,
                   qty=qty1+qty2)
fert2$typ <- ifelse(fert2$typ %in% 1, "DAP", NA)

# -------------------------------------
# read in nitrogen conversion file

conv <- read.csv(file.path(paste0(dataPath,"/../../.."), "Other/Fertilizer/Fert_comp.csv")) %>%
  transmute(typ=Fert_type2, n=N_share/100, p=P_share/100) %>%
  filter(typ %in% c("UREA", "DAP"))

fert1 <- left_join(fert1, conv)
fert2 <- left_join(fert2, conv)

# -------------------------------------
# If purchased amount of nitrogen is zero 
# set to NA to avoid Inf values

fert1$qty <- ifelse(fert1$qty == 0, NA, fert1$qty)
fert2$qty <- ifelse(fert2$qty == 0, NA, fert2$qty)

fert1 <- mutate(fert1,
                Qn=qty*n,
                Qp=qty*p)
fert2 <- mutate(fert2,
                Qn=qty*n,
                Qp=qty*p)

# if Qn is zero change to NA
fert1$Qn <- ifelse(fert1$Qn == 0, NA, fert1$Qn)
fert2$Qn <- ifelse(fert2$Qn == 0, NA, fert2$Qn)

fert <- rbind(fert1, fert2) %>% unique

fert <- group_by(fert, holder_id, household_id, parcel_id, field_id) %>%
  summarise(N=sum(Qn, na.rm=TRUE), P=sum(Qp, na.rm=TRUE),
            UREA=ifelse(any(typ %in% "UREA"), 1, 0),
            DAP=ifelse(any(typ %in% "DAP"), 1, 0))

fert$parcel_id <- as.integer(fert$parcel_id)
fert$field_id <- as.integer(fert$field_id)
fert$WPn <- NA # Empty column to match with 2013
rm(fert1, fert2, conv)

#######################################
############### LABOUR ################
#######################################

# POST HARVEST - unlike 2012-13 survey, there is no
# post planting labour recorded.
ph_lab <- read_dta(file.path(dataPath, "sect10_ph_w1.dta")) %>%
  select(holder_id, household_id, parcel_id, field_id,
                crop_code, ph_s10q01_a:ph_s10q03_f) %>%
  transmute(holder_id, household_id, parcel_id, field_id, crop_code,
            id1=ph_s10q02_a, lab1=ph_s10q02_b*ph_s10q02_c,
            id2=ph_s10q02_e, lab2=ph_s10q02_f*ph_s10q02_g,
            id3=ph_s10q02_i, lab3=ph_s10q02_j*ph_s10q02_k,
            id4=ph_s10q02_m, lab4=ph_s10q02_n*ph_s10q02_o,
            hirM=ph_s10q01_a*ph_s10q01_b,
            hirF=ph_s10q01_d*ph_s10q01_e,
            hirC=ph_s10q01_g*ph_s10q01_h,
            OHHlabM=ph_s10q03_a*ph_s10q03_b,
            OHHlabF=ph_s10q03_c*ph_s10q03_d,
            OHHlabC=ph_s10q03_e*ph_s10q03_f
  )

# -------------------------------------
# make all NA values zero
ph_lab[is.na(ph_lab)] <- 0
ph_lab$crop_code <- as.integer(ph_lab$crop_code)

# sum all labour across a single plot - all measured in days
ph_lab <- transmute(ph_lab, holder_id, household_id, parcel_id, field_id,
                    crop_code, harv_lab=lab1 + lab2 + lab3 + lab4,
                    harv_lab_hire = hirM + hirF + hirC + OHHlabM + OHHlabF + OHHlabC) %>%
         unique

#ph_lab$harv_lab[ph_lab$harv_lab %in% 0] <- NA
#ph_lab$harv_lab_hire[ph_lab$harv_lab_hire %in% 0] <- NA

#######################################
############### GEO ###################
#######################################

geo <- readRDS(file.path(dataPath, "../../../Other/Spatial/ETH/ETH_geo_2011.rds")) 

#######################################
############### AREAs #################
#######################################

# -------------------------------------
# imputed and original gps measurements
# included

areas <- read_dta(paste(dataPath, "../../../Other/Plot_size/areas_eth_y1_imputed.dta", sep="/"))
areas <- select(areas, holder_id, household_id=case_id, parcel_id,
                field_id, area_gps, area_gps_mi_50,
                area_farmer=area_sr)

areas$area_gps <- ifelse(areas$area_gps %in% 0, NA, areas$area_gps)
areas$area_gps_mi_50 <- ifelse(areas$area_gps_mi_50 %in% 0, NA, areas$area_gps_mi_50)

areaTotal <- group_by(areas, household_id) %>%
  summarise(area_tot = sum(area_gps_mi_50, na.rm=TRUE))

areas$parcel_id <- as.integer(areas$parcel_id)

#######################################
########### MISCELLANEOUS #############
#######################################

misc <- remove_all_labels(read_dta(file.path(dataPath, "sect7_pp_w1.dta"))) %>%
  select(household_id, holder_id, rotation = pp_s7q01, credit = pp_s7q06, extension2 = pp_s7q08, oxen = pp_s7q11)

misc$rotation <- ifelse(misc$rotation %in% 2, 0, misc$rotation)
misc$credit <- ifelse(misc$credit %in% 2, 0, misc$credit)
misc$extension2 <- ifelse(misc$extension2 %in% 2, 0, misc$extension2)

#######################################
############## COMMUNITY ##############
#######################################

com3 <- read_dta(file.path(dataPath, "sect3_com_w1.dta")) %>%
  select(ea_id, popEA=cs3q02, HHEA=cs3q03)

com4 <- read_dta(file.path(dataPath, "sect4_com_w1.dta")) %>%
  select(ea_id, road=cs4q01, cost2small_town=cs4q10,
         cost2large_town1=cs4q13_1, cost2large_town2=cs4q13_2,
         bank=cs4q45, micro_finance=cs4q47)

com4$road <- toupper(as_factor(com4$road))
com4$bank <- toupper(as_factor(com4$bank))
com4$micro_finance <- toupper(as_factor(com4$micro_finance))

com6 <- read_dta(file.path(dataPath, "sect6_com_w1.dta")) %>%
  select(ea_id, plant_month1=cs6q03_a, plant_month2=cs6q03_b, plant_month3=cs6q03_c,
         harv_month1=cs6q04_a, harv_month2=cs6q04_b, harv_month3=cs6q04_c,
         ext_agent=cs6q08, dist2ext_agent1=cs6q09_1, dist2ext_agent2=cs6q09_2, fert_source=cs6q12,
         pest_source=cs6q13, seed_source=cs6q14)

# make a seperate file for month data
com6$plant_month1 <- month(com6$plant_month1, label=TRUE)
com6$plant_month2 <- month(com6$plant_month2, label=TRUE)
com6$plant_month3 <- month(com6$plant_month3, label=TRUE)
com6$harv_month1 <- month(com6$harv_month1, label=TRUE)
com6$harv_month2 <- month(com6$harv_month2, label=TRUE)
com6$harv_month3 <- month(com6$harv_month3, label=TRUE)
com6$ext_agent <- as_factor(com6$ext_agent)
com6$fert_source <- as_factor(com6$fert_source)
com6$pest_source <- as_factor(com6$pest_source)
com6$seed_source <- as_factor(com6$seed_source)

com <- left_join(com3, com4) %>% left_join(com6)
rm(com3, com4, com6)

#######################################
########### CROSS SECTION #############
#######################################

# -------------------------------------
# community and location level joins

ETH2011 <- left_join(location, com); rm(com, location)

# -------------------------------------
# household level joins

ETH2011 <- left_join(HH11, ETH2011); rm(HH11) 
ETH2011 <- left_join(ETH2011, areaTotal); rm(areaTotal)

# -------------------------------------
# parcel level joins

ETH2011 <- left_join(ETH2011, parcel); rm(parcel)
ETH2011 <- left_join(ETH2011, misc); rm(misc)
# -------------------------------------
# field level joins

ETH2011 <- left_join(ETH2011, fert); rm(fert) 

# joining with areas causes an attribute error.
# Use Michiel's function to strip attributes from
# areas
stripAttributes <- function(df){
  df[] <- lapply(df, as.vector)
  return(df)
}
areas <- stripAttributes(areas)
ETH2011 <- left_join(ETH2011, areas); rm(areas)
ETH2011 <- left_join(ETH2011, field); rm(field)
ETH2011 <- left_join(ETH2011, geo); rm(geo)

# -------------------------------------
# crop level joins
ETH2011 <- left_join(ETH2011, oput); rm(oput)
ETH2011 <- left_join(ETH2011, crop); rm(crop)
ETH2011 <- left_join(ETH2011, ph_lab); rm(ph_lab)
# attributes also a problem with seed
seed <- stripAttributes(seed)
ETH2011 <- left_join(ETH2011, seed); rm(seed)

# make a surveyyear variable
ETH2011$surveyyear <- 2011

rm(dataPath, make0, x)
ETH2011 <- remove_all_labels(ETH2011)
