###This script imports geographic/shapefile information and WRP characteristics
###Pertaining to the IDPH locations only

library(tidyverse)
library(sf)

################################################################################
#Background Geos (this requires the accompanying files, not just the .shp, to
#be present in the same directory as the .shp file).

#Read in all IL census tracts (from census.gov):
il_tracts <- read_sf("../input/cb_2022_17_tract_500k.shp")

#Join census tracts to create state outline:
IL <- st_union(il_tracts)

################################################################################
#WWTP geographic data (also includes site metadata):
idph_sites <- read_csv("../input/IDPH_sites.csv", 
                       col_names = TRUE,
                       col_types = "cccfffff-fd-ddf--id----f---------") %>%
  filter(!is.na(gps_coordinates_lat) & !is.na(gps_coordinates_long))

#Convert to sf and set CRS:
idph_sf <- st_as_sf(idph_sites, 
                    coords = c("gps_coordinates_long", 
                               "gps_coordinates_lat"), 
                    crs = 4326) %>%
  st_transform(crs = st_crs(il_tracts))

#Add log-transformed variables
idph_sf$log10_pop <- log10(idph_sf$population_served)
idph_sf$log_mgd <- log10(idph_sf$actual_avg_facility_flow_mgd)