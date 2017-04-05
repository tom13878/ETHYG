#'=================================================================================================
#' Project:  IMAGINE ETH
#' Subject:  get data path 
#' Author:   Tom Morley
#' Contact:  tomas.morley@wur.nl
#' Output:   correct datapath for user
#'=================================================================================================

# Use this file to set your path to the data
# used for the Ethiopia LSMS-ISA projects
# check your computer username using
# Sys.info()["user"] and use this in the if
# statement. Then add your dataPath within 
# the {} brackets

# Tom
if(Sys.info()["user"] == "morle001") {
  dataPath <- "C:/Users/morle001/WEcR/ETH"}

# Michiel WEcR
if(Sys.info()["user"] == "dijk158") {
  dataPath <- "C:\\Users\\dijk158\\OneDrive - IIASA\\SurveyData\\ETH"}

# Michiel IIASA
if(Sys.info()["user"] == "") {
  dataPath <- ""}

# Anybody else:
if(Sys.info()["user"] == "") {
  dataPath <- ""}