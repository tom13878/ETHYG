# -------------------------------------
# Ethiopia data 2013 - 2014 survey
# -------------------------------------

if(Sys.info()["user"] == "Tomas"){
  dataPath <- "C:/Users/Tomas/Documents/LEI/data/ETH/2013/Data"
} else {
  dataPath <- "N:/Internationaal Beleid  (IB)/Projecten/2285000066 Africa Maize Yield Gap/SurveyData/ETH/2013/Data"
}

library(haven)
library(dplyr)
library(sjmisc)
library(car)

options(scipen=999)

#######################################
############## LOCATION ###############
#######################################

# file containing the zone, region and district

location <- read_dta(file.path(dataPath, "Household/sect1_hh_w2.dta")) %>%
  select(household_id2, REGCODE = saq01, ZONECODE = saq02, WOREDACODE = saq03, KEBELECODE = saq06, ea_id2, rural) %>%
  unique
location$type <- factor(location$rural, levels=1:3, labels=c("RURAL", "SMALL TOWN", "LARGE TOWN"))
location$rural <- ifelse(location$rural %in% 1, 1, 0)
location$REGCODE <- as.integer(location$REGCODE)

# match up with the names from the survey (prepared in a seperate file)

REGZONE <- read.csv(file.path(paste0(dataPath,"/../../.."), "Other/Spatial/ETH/REGZONEETH.csv"))

# join with household identifications
location <- left_join(location, REGZONE)

rm(REGZONE)

#######################################
########### SOCIO/ECONOMIC ############
#######################################

HH13 <- read_dta(file.path(dataPath, "Household/sect1_hh_w2.dta")) %>%
  select(household_id, household_id2, individual_id, individual_id2,
         ea_id, ea_id2, status=hh_s1q02, sex=hh_s1q03, age=hh_s1q04_a,
         religion=hh_s1q07, marital=hh_s1q08)

HH13$status <- toupper(as_factor(HH13$status))
HH13$religion <- toupper(as_factor(HH13$religion))
HH13$sex <- recode(HH13$sex, c("2=1; 1=0; else = NA"))
HH13$marital <- toupper(as_factor(HH13$marital)) 

# make a new variable cage (cut age = cage) which splits
# individuals according to their age group with
# breaks at 15, 55 and the max age

HH13$cage <- cut(HH13$age, breaks = c(0, 15, 55, max(HH13$age, na.rm=TRUE)),
                 labels=c("0-15", "16-55", "56+"), include.lowest = TRUE, right = TRUE)

# education of household members and sum
# of education of all household members
# between the ages of 15 and 55

education <- read_dta(file.path(dataPath, "Household/sect2_hh_w2.dta")) %>%
  select(household_id2, individual_id, individual_id2, ea_id, ea_id2,
         literate=hh_s2q02, ed_any=hh_s2q03)

education$literate <- ifelse(education$literate %in% 2, 0, education$literate)
education$ed_any <- ifelse(education$ed_any %in% 2, 0, education$ed_any)

# summarise the data: get sum of education
# of household members 15-55 and number of
# household members 15:55

HH13_x <- group_by(HH13, household_id2) %>%
  summarise(N1555=sum(cage %in% "16-55"),
            family_size=n())

# -------------------------------------
# death in the family
# -------------------------------------

death <- read_dta(file.path(dataPath, "Household/sect8_hh_w2.dta")) %>%
  select(household_id2, code=hh_s8q00, death=hh_s8q01) %>% 
  filter(code %in% c("101", "101b")) %>% select(-code) %>%
  group_by(household_id2) %>%
  summarise(death=ifelse(any(death %in% 1), 1, 0))

HH13 <- left_join(HH13, education) %>%
  left_join(HH13_x) %>%
  left_join(death); rm(education, HH13_x, death)

# -------------------------------------
# off-farm income (last 12 months)
# unfortunately we do not see how many
# days are worked per week!!!
# -------------------------------------

off_farm_income <- read_dta(file.path(dataPath, "Household/sect4_hh_w2.dta")) %>%
  select(household_id2, individual_id2, ea_id2,
         main_job=hh_s4q09,
         mj_industry=hh_s4q11_b, mj_employer=hh_s4q12,
         mj_months=hh_s4q13, mj_weekspm=hh_s4q14,
         mj_hourspw=hh_s4q15,mj_wage=hh_s4q16,
         mj_payPeriod=hh_s4q17, mj_grat=hh_s4q18,
         mj_gratPeriod=hh_s4q19,second_job=hh_s4q20,
         sj_industry=hh_s4q22_b, sj_employer=hh_s4q23,
         sj_months=hh_s4q24, sj_weekspm=hh_s4q25,
         sj_hourspw=hh_s4q26, sj_wage=hh_s4q27,
         sj_payPeriod=hh_s4q28, sj_grat=hh_s4q29,
         sj_gratPeriod=hh_s4q30, PSNP=hh_s4q31,
         PSNP_days=hh_s4q32, PSNP_wage=hh_s4q33,
         other_job=hh_s4q34, oj_days=hh_s4q35,
         oj_wage=hh_s4q36)

off_farm_income$mj_pay <- NA
off_farm_income$mj_grat_pay <- NA
off_farm_income$sj_pay <- NA
off_farm_income$sj_grat_pay <- NA
off_farm_income$PSNP_pay <- off_farm_income$PSNP_wage
off_farm_income$oj_pay <- off_farm_income$oj_wage

off_farm_income <- transmute(off_farm_income, household_id2, individual_id2,
                      ea_id2, main_job, mj_industry, mj_employer,
                      mj_pay=ifelse(mj_payPeriod %in% 1, mj_months*mj_weekspm*mj_hourspw*mj_wage, mj_pay),
                      mj_pay=ifelse(mj_payPeriod %in% 3, mj_months*mj_weekspm*mj_wage, mj_pay),
                      mj_pay=ifelse(mj_payPeriod %in% 4, mj_months*mj_weekspm*mj_wage/2, mj_pay),
                      mj_pay=ifelse(mj_payPeriod %in% 5, mj_months*mj_wage, mj_pay),
                      mj_pay=ifelse(mj_payPeriod %in% 6, mj_months*mj_wage/3, mj_pay),
                      mj_pay=ifelse(mj_payPeriod %in% 7, mj_months*mj_wage/6, mj_pay),
                      mj_pay=ifelse(mj_payPeriod %in% 8, mj_months*mj_wage/12, mj_pay),
                      
                      mj_grat_pay=ifelse(mj_gratPeriod %in% 1, mj_months*mj_weekspm*mj_hourspw*mj_grat, mj_grat_pay),
                      mj_grat_pay=ifelse(mj_gratPeriod %in% 3, mj_months*mj_weekspm*mj_grat, mj_grat_pay),
                      mj_grat_pay=ifelse(mj_gratPeriod %in% 4, mj_months*mj_weekspm*mj_grat/2, mj_grat_pay),
                      mj_grat_pay=ifelse(mj_gratPeriod %in% 5, mj_months*mj_grat, mj_grat_pay),
                      mj_grat_pay=ifelse(mj_gratPeriod %in% 6, mj_months*mj_grat/3, mj_grat_pay),
                      mj_grat_pay=ifelse(mj_gratPeriod %in% 7, mj_months*mj_grat/6, mj_grat_pay),
                      mj_grat_pay=ifelse(mj_gratPeriod %in% 8, mj_months*mj_grat/12, mj_grat_pay),
                      
                      second_job, sj_industry, sj_employer,
                      sj_pay=ifelse(sj_payPeriod %in% 1, sj_months*sj_weekspm*sj_hourspw*sj_wage, sj_pay),
                      sj_pay=ifelse(sj_payPeriod %in% 3, sj_months*sj_weekspm*sj_wage, sj_pay),
                      sj_pay=ifelse(sj_payPeriod %in% 4, sj_months*sj_weekspm*sj_wage/2, sj_pay),
                      sj_pay=ifelse(sj_payPeriod %in% 5, sj_months*sj_wage, sj_pay),
                      sj_pay=ifelse(sj_payPeriod %in% 6, sj_months*sj_wage/3, sj_pay),
                      sj_pay=ifelse(sj_payPeriod %in% 7, sj_months*sj_wage/6, sj_pay),
                      sj_pay=ifelse(sj_payPeriod %in% 8, sj_months*sj_wage/12, sj_pay),
                      
                      sj_grat_pay=ifelse(sj_gratPeriod %in% 1, sj_months*sj_weekspm*sj_hourspw*sj_grat, sj_grat_pay),
                      sj_grat_pay=ifelse(sj_gratPeriod %in% 3, sj_months*sj_weekspm*sj_grat, sj_grat_pay),
                      sj_grat_pay=ifelse(sj_gratPeriod %in% 4, sj_months*sj_weekspm*sj_grat/2, sj_grat_pay),
                      sj_grat_pay=ifelse(sj_gratPeriod %in% 5, sj_months*sj_grat, sj_grat_pay),
                      sj_grat_pay=ifelse(sj_gratPeriod %in% 6, sj_months*sj_grat/3, sj_grat_pay),
                      sj_grat_pay=ifelse(sj_gratPeriod %in% 7, sj_months*sj_grat/6, sj_grat_pay),
                      sj_grat_pay=ifelse(sj_gratPeriod %in% 8, sj_months*sj_grat/12, sj_grat_pay),
                      
                      PSNP_wage, oj_wage)

# variable for total pay received from
# main job, secondary job, graduities,
# PSNP and other jobs -> 0s are 
# returned to NAs if necessary

off_farm_income$off_farm_income <- with(off_farm_income,
                     rowSums(cbind(mj_pay,mj_grat_pay,sj_pay, sj_grat_pay, PSNP_wage, oj_wage),
                             na.rm=TRUE))
miss <- with(off_farm_income,
             is.na(mj_pay) & is.na(mj_grat_pay & is.na(sj_pay) & is.na(sj_grat_pay) & is.na(PSNP_wage) & is.na(oj_wage)))
off_farm_income$off_farm_income[miss] <- NA; rm(miss)

# summarise at the household level to get
# a measure of total household off farm
# income

off_farm_income_x <- group_by(off_farm_income, household_id2) %>%
  summarise(off_farm_income_hh=sum(off_farm_income, na.rm=TRUE))

# join off farm income with the household info
HH13 <- left_join(HH13, off_farm_income) %>% left_join(off_farm_income_x)
rm(off_farm_income, off_farm_income_x)
                      
#######################################
############### OUTPUT ################
#######################################

oput <- read_dta(file.path(dataPath, "/Post-Harvest/sect9_ph_w2.dta")) %>%
  select(holder_id, household_id2, parcel_id, field_id,
         crop_code, crop_qty_harv=ph_s9q05, inter_crop=ph_s9q01,
         harv_area=ph_s9q09, harv_month_start=ph_s9q07_a, harv_month_end=ph_s9q07_b)

oput$crop_name <- toupper(as_factor(oput$crop_code))
oput$crop_code <- as.integer(oput$crop_code)
oput$harv_month_start <- toupper(as_factor(oput$harv_month_start))
oput$harv_month_end <- toupper(as_factor(oput$harv_month_end))
oput$inter_crop <- toupper(as_factor(oput$inter_crop))

# -------------------------------------
# add a dummy if a legume was grown
# count number of crops per field
# crop codes taken from appendix in 
# the BID

legumes <- c(11:18, 36, 118)

oput_x <- group_by(oput, holder_id, household_id2, parcel_id, field_id) %>%
  summarise(crop_count=sum(!is.na(crop_code)),
            legume = ifelse(any(crop_name %in% legumes), 1, 0))

oput <- left_join(oput, oput_x); rm(oput_x)

# remove observations with quantity NA or 0
oput <- oput[! is.na(oput$crop_qty_harv) & !oput$crop_qty_harv %in% 0, ]

rm("legumes")

# the price of each crop and how much was
# sold is stored in a seperate section of 
# the household (section 11)

oput2 <- read_dta(file.path(dataPath, "/Post-Harvest/sect11_ph_w2.dta")) %>%
  select(holder_id, household_id2, crop_code,
         sold=ph_s11q01, sold_qty_kg=ph_s11q03_a, sold_qty_gr=ph_s11q03_b,
         value=ph_s11q04, sold_month=ph_s11q06_a, sold_year=ph_s11q06_b,
         trans_cost=ph_s11q09)

oput2$sold_qty_kg[is.na(oput2$sold_qty_k)] <- 0
oput2$sold_qty_gr[is.na(oput2$sold_qty_gr)] <- 0
oput2$sold_qty <- oput2$sold_qty_kg + (oput2$sold_qty_gr/1000)
oput2$sold_qty[oput2$sold_qty==0] <- NA
oput2$sold_qty_kg <-NULL
oput2$sold_qty_gr <-NULL
oput2$crop_price <- oput2$value/oput2$sold_qty
oput2$crop_code <- as.integer(oput2$crop_code)
oput2$sold_month <- toupper(as_factor(oput2$sold_month))
oput2$sold <- toupper(as_factor(oput2$sold))

oput <- left_join(oput, oput2) %>% unique; rm(oput2)

#######################################
########## on farm income #############
#######################################

# from crops
# see output section

# from permanent crops (trees)
tree <- read_dta(file.path(dataPath, "/Post-Harvest/sect12_ph_w2.dta")) %>%
  select(holder_id, household_id2, ea_id2, crop_code,
         sold_tree=ph_s12q06, sold_tree_kg=ph_s12q07, tree_value=ph_s12q08)

# from rented land
rent <- read_dta(file.path(dataPath, "Post-Planting/sect2_pp_w2.dta")) %>%
  select(holder_id, household_id2, ea_id2, parcel_id, rented=pp_s2q10, rented12=pp_s2q11,
         fields_rented=pp_s2q12, rent_cash=pp_s2q13_a, rent_in_kind=pp_s2q13_b)

# from livestock
# information on value of livestock not forthcoming

#######################################
############## CHEMICAL ###############
#######################################

# -------------------------------------
# variables recorded at parcel, field,
# and crop level

# parcel level
parcel <- read_dta(file.path(dataPath, "Post-Planting/sect2_pp_w2.dta")) %>%
  dplyr::select(holder_id, household_id2, parcel_id, fields=pp_s2q02,
                title=pp_s2q04, soil_type=pp_s2q14, soil_qlty=pp_s2q15)

parcel$soil_type <- toupper(as_factor(parcel$soil_type))
parcel$soil_qlty <- toupper(as_factor(parcel$soil_qlty))
parcel$title <- ifelse(parcel$title %in% 2, 0, parcel$title)

# field level variables
# WDswitch

field <- read_dta(file.path(dataPath, "Post-Planting/sect3_pp_w2.dta")) %>%
  dplyr::select(holder_id, household_id2, parcel_id, field_id,
                crop_stand=pp_s3q03b, fallow10=pp_s3q03c, fallow_year=pp_s3q03d,
                extension=pp_s3q11, irrig=pp_s3q12, fert_any=pp_s3q14,
                other_inorg=pp_s3q20a, manure=pp_s3q21, compost=pp_s3q23,
                other_org=pp_s3q25, eros_prot=pp_s3q32, mulch=pp_s3q37) 

field$irrig <- ifelse(field$irrig %in% 2, 0, field$irrig)
field$crop_stand <- toupper(as_factor(field$crop_stand))
field$extension <- ifelse(field$extension %in% 2, 0, field$extension)
field$parcel_id <- as.integer(field$parcel_id)
field$field_id <- as.integer(field$field_id)

# crop level variables
# WDswitch

crop <- read_dta(file.path(dataPath, "/Post-Planting/sect4_pp_w2.dta")) %>%
  dplyr::select(holder_id, household_id2, parcel_id, field_id, crop_code,
                cropping=pp_s4q02, month=pp_s4q12_a, crop_area=pp_s4q03,
                herb=pp_s4q06, fung=pp_s4q07, impr=pp_s4q11, 
                seed_q=pp_s4q11b)

crop$herb <- ifelse(crop$herb %in% 2, 0, crop$herb)
crop$fung <- ifelse(crop$fung %in% 2, 0, crop$fung)
crop$impr <- recode(crop$impr, c("2=1; 1=0; 3=NA"))
crop$crop_code <- as_factor(crop$crop_code)
crop$cropping <- as_factor(crop$cropping)
crop$month <- as_factor(crop$month)
crop$crop_code <- as.integer(crop$crop_code)

# In wave 2 questions about seed are recored at holder level not parcel and field level such as in 2011, hence not useful!
# seed <- remove_all_labels(read_dta(file.path(dataPath, "/Post-Planting/sect5_pp_w2.dta"))) %>%
#   select(household_id, holder_id, parcel_id, field_id, crop_code,
#          impr2 = pp_s5q01, seed_q_kilo = pp_s5q19_a, seed_q_gram = pp_s5q19_b) 
# 
# seed$seed_q_kilo[is.na(seed$seed_q_kilo)] <- 0
# seed$seed_q_gram[is.na(seed$seed_q_gram)] <- 0
# seed$seed_q = seed$seed_q_kilo + seed$seed_q_gram
# seed$seed_q[seed$seed_q == 0] <- NA
# seed$seed_q_kilo <-NULL
# seed$seed_q_gram <-NULL

# -------------------------------------
# unit of observation is not fertilizer

fert1 <- read_dta(file.path(dataPath, "/Post-Planting/sect3_pp_w2.dta")) %>%
  dplyr::select(holder_id, household_id2, parcel_id, field_id, typ=pp_s3q15, qty=pp_s3q16_a,
                purch=pp_s3q16b, purch_kg=pp_s3q16c, valu=pp_s3q16d)
  
fert1$typ <- ifelse(fert1$typ %in% 1, "UREA", NA)
fert1$purch <- ifelse(fert1$purch %in% 1, 1, 0)

fert2 <- read_dta(file.path(dataPath, "/Post-Planting/sect3_pp_w2.dta")) %>%
  dplyr::select(holder_id, household_id2, parcel_id, field_id, typ=pp_s3q18, qty=pp_s3q19_a,
                purch=pp_s3q19b, purch_kg=pp_s3q19c, valu=pp_s3q19d) 

fert2$typ <- ifelse(fert2$typ %in% 1, "DAP", NA)
fert2$purch <- ifelse(fert2$purch %in% 1, 1, 0)

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

fert1$purch_kg <- ifelse(fert1$purch_kg == 0, NA, fert1$purch_kg)
fert2$purch_kg <- ifelse(fert2$purch_kg == 0, NA, fert2$purch_kg)

fert1 <- mutate(fert1,
               Vfert=valu/purch_kg,
               Qn=qty*n,
               Qp=qty*p)
fert2 <- mutate(fert2,
                Vfert=valu/purch_kg,
                Qn=qty*n,
                Qp=qty*p)

# if Qn is zero change to NA
fert1$Qn <- ifelse(fert1$Qn == 0, NA, fert1$Qn)
fert2$Qn <- ifelse(fert2$Qn == 0, NA, fert2$Qn)

# if vfert is 0 change to NA
fert1$Vfert <- ifelse(fert1$Vfert == 0, NA, fert1$Vfert)
fert2$Vfert <- ifelse(fert2$Vfert == 0, NA, fert2$Vfert)

fert <- rbind(fert1, fert2) %>% unique

fert$Pn <- fert$Vfert/fert$n

fert <- group_by(fert, holder_id, household_id2, parcel_id, field_id) %>%
  summarise(N=sum(Qn, na.rm=TRUE), P=sum(Qp, na.rm=TRUE),
            WPn=sum((Qn/N)*Pn, na.rm=TRUE)) %>%
  mutate(WPn = replace(WPn, WPn==0, NA))

rm(fert1, fert2, conv)

#######################################
############### LABOUR ################
#######################################

# POST PLANTING labour
pp_lab <- read_dta(file.path(dataPath, "/Post-Planting/sect3_pp_w2.dta")) %>%
  dplyr::select(holder_id, household_id2, parcel_id, field_id, pp_s3q27_a:pp_s3q29_f) %>%
  transmute(holder_id, household_id2, parcel_id, field_id,
            id1=pp_s3q27_a, lab1=pp_s3q27_b*pp_s3q27_c,
            id2=pp_s3q27_e, lab2=pp_s3q27_f*pp_s3q27_g,
            id3=pp_s3q27_i, lab3=pp_s3q27_j*pp_s3q27_k,
            id4=pp_s3q27_m, lab4=pp_s3q27_n*pp_s3q27_o,
            id5=pp_s3q27_q, lab5=pp_s3q27_r*pp_s3q27_s,
            id6=pp_s3q27_u, lab6=pp_s3q27_v*pp_s3q27_w,
            id7=pp_s3q27_y, lab7=pp_s3q27_z*pp_s3q27_ca,
            hirM=pp_s3q28_a*pp_s3q28_b,
            hirF=pp_s3q28_d*pp_s3q28_e,
            hirC=pp_s3q28_g*pp_s3q28_h,
            OHHlabM=pp_s3q29_a*pp_s3q29_b,
            OHHlabF=pp_s3q29_c*pp_s3q29_d,
            OHHlabC=pp_s3q29_e*pp_s3q29_f
            )

# make all NA values zero
pp_lab[is.na(pp_lab)] <- 0

# sum all labour across a single plot - all measured in days
pp_lab <- transmute(pp_lab, holder_id, household_id2, parcel_id, field_id,
                  plant_lab=lab1 + lab2 + lab3 + lab4 + lab5 + lab6 + lab7,
                  plant_lab_hire = hirM + hirF + hirC + OHHlabM + OHHlabF + OHHlabC)

# presumably if crop was planted then some labour
# was used. Therefore set all 0's for plant_lab
# to NA. UNDONE because labour is split into two components so this does not hold anymore.

#pp_lab$plant_lab[pp_lab$plant_lab %in% 0] <- NA
#pp_lab$plant_lab[pp_lab$plant_lab_hire %in% 0] <- NA

# POST HARVEST
ph_lab <- read_dta(file.path(dataPath, "/Post-Harvest/sect10_ph_w2.dta")) %>%
  dplyr::select(holder_id, household_id2, parcel_id, field_id,
                crop_code, ph_s10q01_a:ph_s10q03_f) %>%
  transmute(holder_id, household_id2, parcel_id, field_id, crop_code,
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

# sum all labour across a single plot - all measured in days
ph_lab <- transmute(ph_lab, holder_id, household_id2, parcel_id, field_id,
                    crop_code, harv_lab=lab1 + lab2 + lab3 + lab4,
                    harv_lab_hire = hirM + hirF + hirC + OHHlabM + OHHlabF + OHHlabC) %>% 
          unique

#ph_lab$harv_lab[ph_lab$harv_lab %in% 0] <- NA
#ph_lab$harv_lab[ph_lab$harv_lab_hire %in% 0] <- NA

#######################################
############### GEO ###################
#######################################

geo <- readRDS(file.path(dataPath, "../../../Other/Spatial/ETH/ETH_geo_2013.rds")) 


#######################################
############### AREAs #################
#######################################

# -------------------------------------
# imputed and original gps measurements
# included

areas <- read_dta(file.path(dataPath, "../../../Other/Plot_size/areas_ETH2013.dta"))
areas <- select(areas, holder_id, household_id2,
                parcel_id, field_id, area_gps, area_gps_mi50,
                area_farmer=area_sr)

areas$area_gps <- ifelse(areas$area_gps %in% 0, NA, areas$area_gps)
areas$area_gps_mi50 <- ifelse(areas$area_gps_mi50 %in% 0, NA, areas$area_gps_mi50)

areaTotal <- group_by(areas, household_id2) %>%
  summarise(area_tot = sum(area_gps_mi50, na.rm=TRUE))

#######################################
########### MISCELLANEOUS #############
#######################################

misc <- remove_all_labels(read_dta(file.path(dataPath, "/Post-Planting/sect7_pp_w2.dta"))) %>%
  select(household_id2, holder_id, rotation = pp_s7q01, credit = pp_s7q06, extension2 = pp_s7q08, oxen = pp_s7q11)

misc$rotation <- ifelse(misc$rotation %in% 2, 0, misc$rotation)
misc$credit <- ifelse(misc$credit %in% 2, 0, misc$credit)
misc$extension2 <- ifelse(misc$extension2 %in% 2, 0, misc$extension2)

#######################################
############## COMMUNITY ##############
#######################################
  
com3 <- read_dta(file.path(dataPath, "Community/sect3_com_w2.dta")) %>%
  select(ea_id2, popEA=cs3q02, HHEA=cs3q03)

com4 <- read_dta(file.path(dataPath, "Community/sect4_com_w2.dta")) %>%
  select(ea_id2, road=cs4q01, cost2small_town=cs4q10,
         cost2large_town=cs4q13, bank=cs4q45, micro_finance=cs4q47)

com4$road <- toupper(as_factor(com4$road))
com4$bank <- toupper(as_factor(com4$bank))
com4$micro_finance <- toupper(as_factor(com4$micro_finance))

com6 <- read_dta(file.path(dataPath, "Community/sect6_com_w2.dta")) %>%
  select(ea_id2, plant_month1=cs6q03_a, plant_month2=cs6q03_b, plant_month3=cs6q03_c,
         harv_month1=cs6q04_a, harv_month2=cs6q04_b, harv_month3=cs6q04_c,
         ext_agent=cs6q08, dist2ext_agent=cs6q09, fert_source=cs6q12,
         pest_source=cs6q13, seed_source=cs6q14)

# make a seperate file for month data
com6$plant_month1 <- as_factor(com6$plant_month1)
com6$plant_month2 <- as_factor(com6$plant_month2)
com6$plant_month3 <- as_factor(com6$plant_month3)
com6$harv_month1 <- as_factor(com6$harv_month1)
com6$harv_month2 <- as_factor(com6$harv_month2)
com6$harv_month3 <- as_factor(com6$harv_month3)
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

ETH2013 <- left_join(location, com); rm(com, location)

# -------------------------------------
# household level joins

ETH2013 <- left_join(HH13, ETH2013); rm(HH13) 
ETH2013 <- left_join(ETH2013, areaTotal); rm(areaTotal)

# -------------------------------------
# parcel level joins

ETH2013 <- left_join(ETH2013, parcel); rm(parcel)
ETH2013 <- left_join(ETH2013, rent); rm(rent)
ETH2013 <- left_join(ETH2013, misc); rm(misc)

# -------------------------------------
# field level joins

ETH2013 <- left_join(ETH2013, fert); rm(fert) 
ETH2013 <- left_join(ETH2013, areas); rm(areas)
ETH2013 <- left_join(ETH2013, pp_lab); rm(pp_lab)
ETH2013 <- left_join(ETH2013, field); rm(field)
ETH2013 <- left_join(ETH2013, geo); rm(geo)

# -------------------------------------
# crop level joins
ETH2013 <- left_join(ETH2013, oput); rm(oput)
ETH2013 <- left_join(ETH2013, crop); rm(crop)
ETH2013 <- left_join(ETH2013, tree); rm(tree)
ETH2013 <- left_join(ETH2013, ph_lab); rm(ph_lab)

# make a surveyyear variable
ETH2013$surveyyear <- 2013

rm(dataPath)
ETH2013 <- remove_all_labels(ETH2013)
