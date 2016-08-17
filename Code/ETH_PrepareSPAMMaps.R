#########################################################################################################
####################################### PROJECT: IPOP ###################################################
#########################################################################################################
# Purpose: Extract data from IFPRI SPAM and sum to regions/district level
#########################################################################################################
#########################################################################################################

# PACKAGES
BasePackages<- c("foreign", "stringr", "gdata", "car", "reshape2", "RColorBrewer", "plyr", "dplyr", "haven", "tidyr")
lapply(BasePackages, library, character.only = TRUE)
SpatialPackages<-c("rgdal", "ggmap", "raster", "rasterVis", "rgeos", "sp", "mapproj", "maptools", "proj4", "gdalUtils")
lapply(SpatialPackages, library, character.only = TRUE)
AdditionalPackages <- c()
lapply(AdditionalPackages, library, character.only = TRUE)


# SET WORKING DIRECTORY
wdPath <- "D:\\Data\\IPOP\\SPAM\\"
setwd(wdPath)

# R SETTINGS
options(scipen=999) # surpress scientific notation
options("stringsAsFactors"=FALSE) # ensures that characterdata that is loaded (e.g. csv) is not turned into factors
options(digits=4)

# SOURCE FILES

# FUNCTIONS
# Obtain country coordinates for target country
Get.country.shapefile.f <- function(iso3c, lev=0, proj = "+proj=longlat +datum=WGS84"){
  
  #download country boundary as spatialpolygonDF (and rewrite as .shp for convenience)
  targetfile <- paste(iso3c, paste("_adm", lev, ".Rdata", sep=""), sep="")
  if(file.exists(paste(countryPath, targetfile, sep="/"))){
    load(paste(countryPath, targetfile, sep="/"))
  } else {
    gadm=getData('GADM', country=iso3c, level=lev, path=countryPath)
  }
  
  # change projection 
  projection <- proj
  country.sp <- spTransform(gadm, CRS(projection))
  return(country.sp)
}

# Function to crop maps to target country
cropmaps.f <-function(file){ 
  setwd("./Raw")
  print(file)
  raster(file) %>%  
    crop(., countryMap) %>%
    writeRaster(., paste(countryPath, paste(iso3c, file, sep="_"), sep="\\"), overwrite=TRUE)
  setwd(wdPath)
}

# SET TARGET COUNTRY
iso3c <-"ETH"

# CREATE COUNTRY FOLDER AND SET WORKING DIRECTORY
countryPath = paste(getwd(), "Processed", iso3c, sep="/") 
if (!file.exists(countryPath)) dir.create(path = countryPath)

# GET COUNTRY MAP
countryMap <- Get.country.shapefile.f(iso3c, lev=1)

# CROP MAIZE MAPS TO TARGET COUNTRY
SPAMfiles <- list.files(paste(wdPath, "RAW", sep="//"))
l_ply(SPAMfiles, cropmaps.f)


# Rename zones using LSMS names
countryMap@data$region_lsms <- countryMap@data$NAME_1
countryMap@data$region_lsms[countryMap@data$NAME_1 %in% c("Oromia")] <- "Oromiya"
countryMap@data$region_lsms[countryMap@data$NAME_1 %in% c("Somali")] <- "Somalie"
countryMap@data$region_lsms[countryMap@data$NAME_1 %in% c("Benshangul-Gumaz")] <- "Benishangul Gumuz"
countryMap@data$region_lsms[countryMap@data$NAME_1 %in% c("Southern Nations, Nationalities and Peoples")] <- "SNNP"
countryMap@data$region_lsms[countryMap@data$NAME_1 %in% c("Gambela Peoples")] <- "Gambella"
countryMap@data$region_lsms[countryMap@data$NAME_1 %in% c("Harari People")] <- "Harari"
countryMap@data$region_lsms <- factor(countryMap@data$region_lsms)

# New zone names
p.df <- data.frame(ID=1:length(countryMap), zone =  countryMap@data$region_lsms)

# Plot map
colors=rainbow(length(countryMap$region_lsms))
plot(countryMap,col=colors)
#legend("topleft",legend=countryMap$region_lsms,fill=colors,cex=1.3,bty="n" )
text(coordinates(countryMap), labels = countryMap$region_lsms)

# Calculate total (irrigated plus rainfed) production per zone
prodZone <-raster(paste(".\\Processed\\", iso3c, "\\ETH_spam2005v2r0_production_maize_total.tif", sep="")) %>%
                  raster::extract(., countryMap, df=T) %>% rename(Production = ETH_spam2005v2r0_production_maize_total) %>%
                  group_by(ID) %>%
                  summarize(Production = sum(Production, na.rm=T)) %>%
                  left_join(., p.df)

# Calculate weighted total yield per zone
Zone_ETH <-stack(paste(".\\Processed\\", iso3c, "\\ETH_spam2005v2r0_harvested-area_maize_total.tif", sep=""), 
                  paste(".\\Processed\\", iso3c, "\\ETH_spam2005v2r0_yield_maize_total.tif", sep=""),
                  paste(".\\Processed\\", iso3c, "\\ETH_spam2005v2r0_production_maize_total.tif", sep="")) %>%
  raster::extract(., countryMap, df=T) %>%
  dplyr::rename(area = ETH_spam2005v2r0_harvested.area_maize_total, 
                yield = ETH_spam2005v2r0_yield_maize_total,
                production = ETH_spam2005v2r0_production_maize_total) %>%
  group_by(ID) %>%
  summarize(production = sum(production, na.rm=T), yield = sum(area *yield, na.rm=T)/sum(area, na.rm=T)) %>%
  left_join(., p.df)

# Calculate total maize production
# Note that SPAM is calibrated on FAO data for 2004-2006, not 2005.
prodTotalSPAM <- sum(Zone_ETH$production, na.rm=T)

# Load maize production data from FAOSTAT 
# Better to download one file for all countries
prodFAO <- read.csv(paste("Processed", iso3c, "ETH_MaizeProduction.csv", sep="//"))

# Compute scaling factor between SPAM and TZA survey period: 2013 survey takes place in these years (2014 not available in FAOSTAT)
scalingFactor <- prodFAO %>% 
              filter(Year %in% c(2013)) %>%
              summarize(scalingFactor = mean(Value, na.rm=T)/prodTotalSPAM)

scalingFactor <- scalingFactor$scalingFactor

# Compute total production per zone for target period
Zone_ETH$TargetProduction <- Zone_ETH$production*scalingFactor

# Write file
saveRDS(Zone_ETH, file= "D:\\Data\\Projects\\ETHYG\\Cache\\SPAMdata_ETH.rds")
 