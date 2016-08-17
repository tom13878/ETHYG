# PROJECT: IPOP
# ``````````````````````````````````````````````````````````````````````````````````````````````````
# ``````````````````````````````````````````````````````````````````````````````````````````````````
# # http://www.isric.org/data/AfSoilGrids250m
# Purpose: Download spatial precipitation, soil, etc data and prepare country files for yield analysis
# ``````````````````````````````````````````````````````````````````````````````````````````````````
# `````````````````````````````````````````````````````````````````````````````````````````````````` 

# PACKAGES
BasePackages<- c("foreign", "stringr", "gdata", "car", "reshape2", "RColorBrewer", "plyr", "dplyr")
lapply(BasePackages, library, character.only = TRUE)
SpatialPackages<-c("rgdal", "ggmap", "raster", "rasterVis", "rgeos", "sp", "mapproj", "maptools", "proj4", "gdalUtils")
lapply(SpatialPackages, library, character.only = TRUE)
AdditionalPackages <- c("plotKML", "RCurl", "XML", "R.utils", "GSIF")
lapply(AdditionalPackages, library, character.only = TRUE)

#'NOTE THAT
#'The code below requires GDAL to be installed on your computer and probably
#'also SAGA and some other software. The following will identify the paths
#'of the additional software and indicate availability.
#'
# Set working environment
plotKML.env(silent = FALSE)
gdal_setInstallation()

# SET WORKING DIRECTORY
wdPath<-"D:\\Data\\IPOP\\AFSIS\\"
setwd(wdPath)
 
# # # DOWNLOAD AFRICA SOILGRIDS
# # # Download Soilgrids: http://gsif.isric.org/doku.php?id=wiki%3Atutorial_soilgrids
# # # soilgrids:soilgrids are password and username
# sg.ftp <- "ftp://soilgrids:soilgrids@ftp.soilgrids.org/data/AF/recent/"
# filenames = getURL(sg.ftp, ftp.use.epsv = FALSE, dirlistonly = TRUE, userpwd = "soilgrids")
# filenames = strsplit(filenames, "\r*\n")[[1]]
# targetnames <- c("ORCDRC", "BLD", "CRFVOL", "PHIHOX")
# targetdepth <- c("sd1", "sd2", "sd3", "sd4", "sd5", "sd6")
# targetfiles <-filenames[grep(pattern=paste(targetnames, collapse="|"), filenames)]
# targetfiles <-targetfiles[grep(pattern=paste(targetdepth, collapse="|"), targetfiles)]
# l_ply(targetfiles,function(x) if(!file.exists(x)){try(download.file(paste(sg.ftp, x, sep=""),x))})
# 
# # # # unzip tif files
# tifFiles<-targetfiles[grep(pattern=".gz", targetfiles)]
# l_ply(tifFiles, gunzip, skip=TRUE)
# tifFiles<-llply(tifFiles,function(x) strsplit(x,".gz")[[1]][1])
# llply(tifFiles,GDALinfo)

# Download 
# SUBSET FILES FOR A TARGET COUNTRY AND WRITE AS RASTERSTACK
# Set target country
iso3c <-"ETH"

# CREATE COUNTRY FOLDER AND SET WORKING DIRECTORY
countryPath = paste(getwd(), "Processed", iso3c, sep="/") 
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

# Extract information for target country [GDALWARP does not work so different approach applied below]
# AfricaSoilFiles <- list.files(path = wdpath, pattern = "\\.tif$") 
# CountrySoilData <- lapply(AfricaSoilFiles, function(x) {
#                   gdalwarp(srcfile=paste(wdpath, x, sep="\\"), dstfile=paste(iso3c, x, sep="_"), 
#                   t_srs="+proj=longlat +datum=WGS84",
#                   te=coord, overwrite=TRUE, output_Raster=TRUE, verbose=TRUE)
#                 }
#               ) %>% do.call(stack,.)

# # Write raster as brick
# writeRaster(CountrySoilData, paste(countrypath, paste(iso3c, "_SoilData", sep=""), sep="\\", overwrite=TRUE)
# CountrySoilData <- brick(paste(countrypath, paste(iso3c, "_SoilData", sep=""), sep="\\"))

# Reproject country map
country.map <- Get.country.shapefile.f(iso3c)
ProjectionAFSIS <- "+proj=laea +lat_0=5 +lon_0=20 +x_0=0 +y_0=0 +units=m +ellps=WGS84 +datum=WGS84"
Baseprojection <- "+proj=longlat +datum=WGS84"
country_AFSIS_projection <- spTransform(country.map, CRS(ProjectionAFSIS))

# Crop and reproject country files
#' NOTE that gdal and raster that are applied in this code will
#' create large temporary files (several GB) on you computer. In my case:
#' C:\Users\dijk158\AppData\Local\Temp\ [RtmpKsUuSY for GDAL and another for Raster]
#' Make sure there is enough space and if needed delete afterwards.

# Function to crop 250m maps
soilcropreproject.f <-function(file){ 
  setwd("./Raw")
  print(file)
  raster(file) %>%  
    crop(., country_AFSIS_projection) %>% # Perhaps aggregate to 1 km to 1km here
    reproject(., program = "GDAL") %>%
    writeRaster(., paste(countryPath, paste(iso3c, file, sep="_"), sep="\\"), overwrite=TRUE)
  setwd(wdPath)
}

# Function to aggregate 250m maps to 1km and crop
soilcropreproject2.f <-function(file){ 
  setwd("./Raw")
  print(file)
  newfile <- str_sub(file, 0, -9)
  raster(file) %>%  
    crop(., country_AFSIS_projection) %>%
    aggregate(., 4, fun="mean") %>%
    reproject(., program = "GDAL") %>%
    writeRaster(., paste(countryPath, paste(iso3c, newfile, "1km.tif", sep="_"), sep="\\"), overwrite=TRUE)
  setwd(wdPath)
}

soilFiles250 <- list.files(paste(wdPath, "RAW", sep="//"), pattern="_250m.tif$")
l_ply(soilFiles250, soilcropreproject2.f)

# Root depth map
soilcropreproject.f("af_agg_ERZD_TAWCpF23mm__M_1km.tif")

# RootDepth <- raster("D:\\Data\\IPOP\\AFSIS\\Processed\\ETH\\ETH_af_agg_ERZD_TAWCpF23mm__M_1km.tif")
# plot(RootDepth)
# plot(country.map, add=T)
