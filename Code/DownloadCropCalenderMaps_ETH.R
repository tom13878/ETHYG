# PROJECT: DFID
# ``````````````````````````````````````````````````````````````````````````````````````````````````
# ``````````````````````````````````````````````````````````````````````````````````````````````````
# R code to prepare crop calendar data
# ``````````````````````````````````````````````````````````````````````````````````````````````````
# `````````````````````````````````````````````````````````````````````````````````````````````````` 

# INSTALL PACKAGES AND SET WORKING DIRECTORY
BasePackages <- c("foreign", "stringr", "gdata", "car", "reshape2", "RColorBrewer", "plyr", "dplyr")
lapply(BasePackages, library, character.only = TRUE)
SpatialPackages <- c("rgdal", "gdalUtils", "ggmap", "raster", "rasterVis", "rgeos", "sp", "mapproj", "maptools", "proj4")
lapply(SpatialPackages, library, character.only = TRUE)
AdditionalPackages <- c("ncdf")
lapply(AdditionalPackages, library, character.only = TRUE)

# OPTIONS
options(scipen = 999) # surpress scientific notation
options("stringsAsFactors" = FALSE) 
options(digits = 4)

# SET WORKING DIRECTORY
# Note that rainfall data will be saved automatically in the wd.
wdPath <- "D:\\Data\\IPOP\\CropCalendar"
setwd(wdPath)

# SOURCE FUNCTIONS

# SUBSET FILES FOR A TARGET COUNTRY AND WRITE AS RASTERSTACK
# Set target country
iso3c <-"ETH"

# CREATE COUNTRY FOLDER AND SET WORKING DIRECTORY
countryPath = paste(getwd(), "/", iso3c, sep="") 
if (!file.exists(countryPath)) dir.create(path = countryPath)

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

country.map <- Get.country.shapefile.f(iso3c, 1)


# Get crop calendar data
# Get general information from nc file
check <- open.ncdf(".\\Raw\\Maize.crop.calendar.fill.nc")
print(check)
close.ncdf(check)

get_nc_f <- function(nc, var, proj = "+proj=longlat +datum=WGS84"){
  rasterfile <- raster(nc, varname = var)
  projection(rasterfile) <- proj
  return(rasterfile)
}

# Main growing season
plant_start <- get_nc_f(".\\Raw\\Maize.crop.calendar.fill.nc", "plant.start")
plant_end <- get_nc_f(".\\Raw\\Maize.crop.calendar.fill.nc", "plant.end")
harvest_start <- get_nc_f(".\\Raw\\Maize.crop.calendar.fill.nc", "harvest.start")
harvest_end <- get_nc_f(".\\Raw\\Maize.crop.calendar.fill.nc", "harvest.end")

cc <- stack(plant_start, plant_end, harvest_start, harvest_end)
cc
plot(cc)
plot(cc[[1]])
map(add=T)
writeRaster(cc, ".\\Processed\\CropCalendar", overwrite=TRUE)

# Second growing season
plant_start2 <- get_nc_f(".\\Raw\\Maize.2.crop.calendar.fill.nc", "plant.start")
plant_end2 <- get_nc_f(".\\Raw\\Maize.2.crop.calendar.fill.nc", "plant.end")
harvest_start2 <- get_nc_f(".\\Raw\\Maize.2.crop.calendar.fill.nc", "harvest.start")
harvest_end2 <- get_nc_f(".\\Raw\\Maize.2.crop.calendar.fill.nc", "harvest.end")

cc2 <- stack(plant_start2, plant_end2, harvest_start2, harvest_end2)
cc2
plot(cc2)
plot(cc2[[1]])
map(add=T)
writeRaster(cc, ".\\Processed\\CropCalendar", overwrite=TRUE)



# Crop to country
cc.country <- crop(cc, country.map)
plot(cc.country[[3]])
map(add=T)
plot(country.map, add=T)

cc2.country <- crop(cc2, country.map)
plot(cc2.country[[3]])
#map(add=T)
plot(country.map, add=T)
