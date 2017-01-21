#'========================================================================================================================================
#' Project:  IMAGINE
#' Subject:  Script to create maps
#' Author:   Michiel van Dijk
#' Contact:  michiel.vandijk@wur.nl
#'========================================================================================================================================

### PACKAGES
if(!require(pacman)) install.packages("pacman")
# Key packages
p_load("tidyverse", "readxl", "stringr", "scales", "RColorBrewer", "rprojroot")
# Spatial packages
p_load("rgdal", "ggmap", "raster", "rasterVis", "rgeos", "sp", "mapproj", "maptools", "proj4", "gdalUtils")
# Additional packages
#p_load("WDI", "countrycode")

### DETERMINE ROOT PATH
root <- find_root(is_rstudio_project)

### DATAPATH
dataPath <- "C:\\Users\\dijk158\\OneDrive - IIASA\\SurveyData"
GYGApath <- "D:\\Data\\IPOP\\GYGA\\"

### R SETTINGS
options(scipen=999) # surpress scientific notation
options("stringsAsFactors"=FALSE) # ensures that characterdata that is loaded (e.g. csv) is not turned into factors
options(digits=4)


### GYGA MAP PREPARATION

# Load polygon map and cut out tartet country
GYGApath <- "D:\\Data\\IPOP\\GYGA\\"
dsn=paste(GYGApath, "\\CZ_SubSaharanAfrica\\CZ_AFRGYGACNTRY.shp", sep="")
ogrListLayers(dsn)
ogrInfo(dsn, layer="CZ_AFRGYGACNTRY")
GYGA.Africa <- readOGR(dsn, layer = "CZ_AFRGYGACNTRY")
projection(GYGA.Africa) # check projection
GYGA.Africa <- spTransform(GYGA.Africa, CRS("+proj=longlat +datum=WGS84")) 

# Cut out GHA from GYGA map
GYGA_country <- GYGA.Africa[GYGA.Africa$REG_NAME=="Ethiopia",]
plot(GYGA_country)
rm(GYGA.Africa)

# Get GYGA data and add YG
GYGA_yg_df <- read_excel(file.path(GYGApath, "GygaEthiopia.xlsx"), sheet=3) %>%
  dplyr::filter(CROP=="Rainfed maize") %>%
  mutate(YG = (YW-YA)/YW*100,
         CLIMATEZONE = as.character(CLIMATEZONE))

# Link yield gap data
# in order to merge data with spatialpolygondataframe the row.names need to be the same.
# For this reason I first merge the additionald data and then create a spatialpolygondataframe
# again ensuring that the rownames are preserved.

GYGA_df <- as(GYGA_country, "data.frame") %>%
  rename(CLIMATEZONE = GRIDCODE)
GYGA_df$id <-row.names(GYGA_df)
GYGA_df <- left_join(GYGA_df, GYGA_yg_df)
row.names(GYGA_df) <- GYGA_df$id
GYGA_country <- SpatialPolygonsDataFrame(as(GYGA_country, "SpatialPolygons"),data=GYGA_df)

# Maps with GYGA yield potential and plot information
# transform shapefile in dataframe for ggplot. rownames are used as ID for the countries. Using ID gives strange results. 
# Additional data is linked back again
GYGA_df_f <- fortify(GYGA_country) %>%
  left_join(., GYGA_df)

### GYGA YIELD GAP 
# Draw map
Map_GYGA_YG <- ggplot()+
  geom_polygon(data=GYGA_df_f, aes(x=long, y=lat, group=group, fill=YG), colour="black")+
  geom_polygon(data=subset(GYGA_df_f, is.na(YG)), aes(x=long, y=lat, group=group), fill="white", colour="black") +
  #scale_fill_gradientn(colours = terrain.colors(10)) +
  #scale_fill_gradient(low = "light green", high = "dark green") +
  scale_fill_distiller(palette = "Spectral", name = "%") +
  coord_equal() +
  labs(
    title = "Water-limited yield gap in Ethiopia (%)",
    #subtitle = "check",
    caption = "Source: Global Yield Gap Atlas (www.yieldgap.org)",
    x="", y="") +
  theme_classic() +
  theme(legend.key=element_blank(),
        line = element_blank(),
        axis.text = element_blank())
#GYGA_YG



