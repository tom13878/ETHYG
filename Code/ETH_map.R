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


### LOAD DATA
db9 <- readRDS(file.path(root, "Cache/db9.rds"))
db_sfaCD_CRE_Z <- readRDS(file.path(root, "Cache/db_sfaCD_CRE_Z.rds"))
db1 <- readRDS(file.path(root, "Cache/db1.rds"))

### LOAD GYGA MAP
dsn=paste(GYGApath, "\\CZ_SubSaharanAfrica\\CZ_AFRGYGACNTRY.shp", sep="")
ogrListLayers(dsn)
ogrInfo(dsn, layer="CZ_AFRGYGACNTRY")
GYGA.Africa<-readOGR(dsn, layer = "CZ_AFRGYGACNTRY")
projection(GYGA.Africa) # check projection
GYGA.Africa <- spTransform(GYGA.Africa, CRS("+proj=longlat +datum=WGS84"))

# Get GYGA
GYGA.country.yield.data <- read_excel(paste(GYGApath, "GygaEthiopia.xlsx", sep="\\"), sheet=3)

# Select data for maize
GYGA.country.yield.data <- filter(GYGA.country.yield.data, CROP=="Rainfed maize")

# Cut out ETH from GYGA map
GYGA.country <- GYGA.Africa[GYGA.Africa$REG_NAME=="Ethiopia",]

# Link yield gap data
# in order to merge data with spatialpolygondataframe the row.names need to be the same.
# For this reason I first merge the additionald data and then create a spatialpolygondataframe
# again ensuring that the rownames are preserved.

GYGA.country.data <- as(GYGA.country, "data.frame")
GYGA.country.data$id <-row.names(GYGA.country.data)
GYGA.country.data <- merge(GYGA.country.data, GYGA.country.yield.data[,c(1:8)], by.x=c("GRIDCODE"), by.y=c("CLIMATEZONE"), all.x=TRUE, sort=FALSE)
row.names(GYGA.country.data) <- GYGA.country.data$id
GYGA.country <- SpatialPolygonsDataFrame(as(GYGA.country, "SpatialPolygons"),
                                         data=GYGA.country.data)

# Maps with GYGA yield potential and plot information
# transform shapefile in dataframe for ggplot. rownames are used as ID for the countries. Using ID gives strange results. 
# Additional data is linked back again
GYGA.country.fort<- fortify(GYGA.country) 
GYGA.country.fort <- merge(GYGA.country.fort, GYGA.country.data, by="id")
GYGA.country.fort$yieldclass <- cut(GYGA.country.fort$YW, breaks=c(6, 8.5, 11, 13.5, 16, 19))

# GYGA map
Map_GYGA <- ggplot()+
  geom_polygon(data = filter(GYGA.country.fort, !is.na(YW)), aes(x = long, y = lat, group = group, fill = yieldclass), colour="black")+
  geom_polygon(data = filter(GYGA.country.fort, is.na(YW)), aes(x = long, y = lat, group = group), fill="white", colour="black") +
  scale_fill_discrete(name = "Potential water\nlimited yield (tons)") +
  coord_equal()+
  labs(x="", y="")+
  theme_classic() +
  theme(legend.key=element_blank(),
        line = element_blank(),
        axis.text = element_blank())

Map_GYGA

# Prepare mean yield data
yld <- db9 %>%
  ungroup() %>%
  group_by(lon, lat) %>%
  summarize(n=n(),
            av_yld = (sum(Y*area)/sum(area))/1000) %>%
  mutate(av_yld2 = cut(av_yld, breaks=c(0, 1, 2, 3, 10)))

  
# Combined GYGA_LSMS map
Map_GYGA_LSMS <- Map_GYGA +
  geom_point(data = yld, aes(x = lon, y = lat, size = av_yld2), colour = "black")+
  scale_size_manual(name =  "Average yield (tons)", values = c(1, 2, 3, 4)) 

Map_GYGA_LSMS 

# # Zonal map with community yield levels
# # read in map of Ethiopia as a SpatialPolygonsDataFrame
# 
# countryMap <- getData('GADM', country = "ETH", level = 1) 
# 
# # Rename zones using LSMS names
# countryMap@data$ZONE <- countryMap@data$NAME_1
# countryMap@data$ZONE[countryMap@data$NAME_1 %in% c("Oromia")] <- "Oromiya"
# countryMap@data$ZONE[countryMap@data$NAME_1 %in% c("Somali")] <- "Somalie"
# countryMap@data$ZONE[countryMap@data$NAME_1 %in% c("Benshangul-Gumaz")] <- "Benishangul Gumuz"
# countryMap@data$ZONE[countryMap@data$NAME_1 %in% c("Southern Nations, Nationalities and Peoples")] <- "SNNP"
# countryMap@data$ZONE[countryMap@data$NAME_1 %in% c("Gambela Peoples")] <- "Gambella"
# countryMap@data$ZONE[countryMap@data$NAME_1 %in% c("Harari People")] <- "Harari"
# countryMap@data$ZONE <- factor(countryMap@data$ZONE)
# 
# # Remove Addis Ababa and Dire Dawa
# countryMap <- countryMap[!(countryMap@data$ZONE %in% c("Addis Abeba", "Dire Dawa")),]
# plot(countryMap)
# 
# #  fortify spatial data to use with ggplot and join using join functions from dplyr
# #    The join is on id, make sure all ids are character vectors
# tf <- fortify(countryMap)
# countryMap@data <- rename(countryMap@data, id = ID_1)
# countryMap@data$id <- as.character(countryMap@data$id)
# tf2 <- left_join(tf, countryMap@data)
# 
# # Use ggplot to plot map of Tanzania, specifying the labels and choosing nice colours
# #    from the RColorBrewer package
# 
# #display.brewer.all()
# 
# ZONE_LSMS <- ggplot()+
#   geom_polygon(data=tf2, aes(x=long, y=lat, group=group, fill=ZONE), colour="black")+
#   geom_point(data=meanYield, aes(x=lon, y=lat, size=(meanYield2)), colour="black")+
#   scale_fill_brewer(name = "Zones", palette = "Set1") +
#   scale_size_manual(name="Average yield (tons)", values=c(1.5, 2.5, 3.5, 4,5)) +
#   coord_equal()+
#   labs(x="", y="")+
#   theme_classic() +
#   theme(legend.key=element_blank(),
#         line = element_blank(),
#         axis.text = element_blank())
# #ZONE_LSMS 
# 
# 
# 
