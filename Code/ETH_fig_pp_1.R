#'========================================================================================================================================
#' Project:  IMAGINE
#' Subject:  Script to create figures
#' Author:   Michiel van Dijk
#' Contact:  michiel.vandijk@wur.nl
#'========================================================================================================================================

### PACKAGES
if(!require(pacman)) install.packages("pacman")
# Key packages
p_load("tidyverse", "readxl", "stringr", "scales", "RColorBrewer", "rprojroot")
# Spatial packages
#p_load("rgdal", "ggmap", "raster", "rasterVis", "rgeos", "sp", "mapproj", "maptools", "proj4", "gdalUtils")
# Additional packages
p_load("cowplot")

### DETERMINE ROOT PATH
root <- find_root(is_rstudio_project)

### DATAPATH

### R SETTINGS
options(scipen=999) # surpress scientific notation
options("stringsAsFactors"=FALSE) # ensures that characterdata that is loaded (e.g. csv) is not turned into factors
options(digits=4)

### SELF SUFFICIENCY FIGURE
# Load data
ss_ratio_raw <- read_excel(file.path(root, "Data/PNAS_ETH_processed.xlsx"), sheet = "ss_ratio") %>%
  mutate(scenario = factor(scenario))

# Graph
col = rev(brewer.pal(9,"PuBuGn")[c(3,5,7,9)])
ss2010 <- round(ss_ratio_raw$ss_ratio[ss_ratio_raw$type == "Line"], 2)

p_ss = ggplot() +
  theme_classic(16) +
  geom_bar(data = filter(ss_ratio_raw, type == "Bar"), aes(x = type, y = ss_ratio, fill = scenario), colour = "black", stat = "identity", position = "dodge", alpha = 1) +
  scale_y_continuous(limits = c(0,2), expand=c(0,0), breaks = c(seq(0, 2, 0.5), 1, ss2010)) +
  scale_fill_manual(values = col) +
  labs(x = "", y = "Self-sufficiency ratio") +
  theme(axis.title.x=element_blank(),
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      axis.text.y=element_text(colour="black")) +
  geom_hline(aes(yintercept = 1), linetype = "dashed") +
  geom_hline(aes(yintercept = ss2010), linetype = "dashed") +
  annotate("text", 0.55, ss2010+0.02, label = "Self-sufficiency ratio 2010", hjust = 0, size = 4) +
  annotate("text", 0.55, 1.02, label = "Self-suffient", hjust = 0, size = 4) 


### AREA FIGURE
# Load data
ss_area_raw <- read_excel(file.path(root, "Data/PNAS_ETH_processed.xlsx"), sheet = "ss_area") %>%
  na.omit()

# Graph

# Graph
#col = brewer.pal(9,"PuBuGn")[c(3,5,7,9)]
a2010 <- round(ss_area_raw$area[ss_area_raw$scenario == "Current cereal area (Mha)"], 2)
aPotential <- round(ss_area_raw$area[ss_area_raw$scenario == "Potential available cereal area (Mha)"], 2)

p_area = ggplot() +
  theme_classic(16) +
  geom_bar(data = filter(ss_area_raw, type == "Bar"), aes(x = type, y = area, fill = scenario), colour = "black", stat = "identity", position = "dodge", alpha = 1) +
  scale_y_continuous(limits = c(0, 15), expand=c(0,0), breaks = c(seq(0, 15, 5), a2010, aPotential)) +
  scale_fill_manual(values = col) +
  labs(x = "", y = "Cereal area (Mha)") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_text(colour="black")) +
  geom_hline(aes(yintercept = a2010), linetype = "dashed") +
  geom_hline(aes(yintercept = aPotential), linetype = "dashed") +
    annotate("text", 1.45, aPotential+0.25, label = "Potential available cereal area", hjust = 1, size = 4) +
    annotate("text", 1.45, a2010+0.25, label = "Current cereal area", hjust = 1, size = 4)  

### COMBINE GRAPHS
p_tot <- plot_grid(p_ss + theme(legend.position="none"), 
                   p_area + theme(legend.position="none"), labels = c("A", "B"), align = "h")
legend <- get_legend(p_area + theme(legend.position="bottom") + 
                        guides(fill= guide_legend("", ncol=2)))

Fig_ss_a <- plot_grid(p_tot, legend, ncol = 1, rel_heights = c(1, .2))
Fig_ss_a






