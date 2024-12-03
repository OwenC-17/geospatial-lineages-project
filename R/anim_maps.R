library(tidyverse)
library(sf)
library(raster)
library(dplyr)
library(spData)
library(spDataLarge)
library(tmap)
library(leaflet)
library(mapview)
library(ggplot2)
library(shiny)
library(terra)
library(zoo)
library(viridis)
library(naniar)
library(simputation)

mwrd_map <- st_read("../input/MWRD_Files/MWRD.shp")

#Confirm:
#plot(mwrd_map)
#view(mwrd_map)

chicago_sites <- read.csv("../input/chicago_sites.csv", header = FALSE, col.names = c("number", "sites")) %>%
  mutate(sites = gsub('\\"', '', sites))

long_ww_chicago <- filter(long_ww_lin_w_sample_info, 
                          wwtp_name %in% chicago_sites$sites) %>%
  mutate(wwtp_name = gsub('METROPOLITAN ', 'M', wwtp_name))


long_ww_chicago$NAME <- recode(as.factor(long_ww_chicago$wwtp_name),
                               "MWRDGC-KIRIE" = "Kirie",
                               "MWRDGC HANOVER PARK WRP" = "Hanover",
                               "MWRDGC STICKNEY WRP - SWS(N/S)" = "Stickney",
                               "MWRDGC TERRENCE J O'BRIEN WTR RECLAMATION PLANT" = "OBrien",
                               "MWRDGC-JOHN E. EGAN WRP" = "Egan",
                               "MWRDGC LEMONT WRP" = "Lemont",
                               "MWRDGC CALUMET WRP" = "Calumet",
                               "MWRDGC STICKNEY WRP - WS" = "Stickney")

mwrd_map_joined <- left_join(mwrd_map, long_ww_chicago, by = "NAME")

ba1_map <- mwrd_map_joined %>%
  filter(named_variant_id == "Omicron BA.1")

ba1_anim = tm_shape(ba1_map) + tm_polygons(col = "abvec") +
  tm_facets(along = "weeks_since_start", free.coords = FALSE) #using "along =" instead of "by =" generates animation. "free.coords = FALSE" keeps the viewing window stationary.

ba5_map <- mwrd_map_joined %>%
  filter(named_variant_id == "Omicron BA.5") %>%
  group_by(wwtp_name) %>%
  mutate(abvec = ifelse(is.na(abvec), 0, abvec)) %>%
  mutate(rollav = rollmean(abvec, k = 5, fill = 0, na.pad = TRUE, align = 'right'))

ba5_map_grouped <- ba5_map %>%
  group_by(wwtp_name, weeks_since_start, sample_collect_date, named_variant_id) %>%
  summarize(totalBa5 = sum(abvec))

ba5_anim = tm_shape(ba5_map) + tm_borders + tm_polygons(col = "rollav", palette = "viridis") +
  tm_facets(along = "weeks_since_start", free.coords = FALSE) #using "along =" instead of "by =" generates animation. "free.coords = FALSE" keeps the viewing window stationary.


ba5_anim2 = tm_shape(ba5_map_grouped) + tm_borders + tm_polygons(col = "totalBa5") +
  tm_facets(along = "weeks_since_start", free.coords = FALSE) #using "along =" instead of "by =" generates animation. "free.coords = FALSE" keeps the viewing window stationary.


ba5_anim2

tmap_animation(ba5_anim, filename = "ba5_anim.gif", delay = 25)




####Miscellaneous exploratory stuff
table(mwrd_map_joined$named_variant_id)


long_ww_chicago %>%
  filter(NAME == "Kirie") %>%
  ggplot(aes(x = sample_collect_date, y = abvec)) + geom_area(position = "fill", aes(fill = named_variant_id))


head(mwrd_map_joined)


#Interpolation
###############################################################################
combos <- mwrd_map_joined %>%
  complete(wwtp_name, sample_collect_date, named_variant_id)
grouped <- combos %>%
  group_by(wwtp_name, sample_collect_date, named_variant_id)

voc_sums <- grouped %>%
  summarize(voc_abundance = sum(abvec)) %>%
  mutate(named_variant_id = factor(named_variant_id)) %>%
  mutate(wwtp_name = factor(wwtp_name)) %>%
  filter(voc_abundance <= 1 | is.na(voc_abundance))

filled <- voc_sums %>%
  bind_shadow() %>%
  mutate(numdat = as.numeric(sample_collect_date)) %>%
  group_by(named_variant_id, wwtp_name) %>%
  fill(voc_abundance, .direction = "down")

ggplot(filled, aes(x = sample_collect_date, y = voc_abundance)) + 
  geom_point(aes(col = named_variant_id, shape = voc_abundance_NA)) + facet_wrap(~wwtp_name)

ggplot(filter(filled, named_variant_id == "Delta" & wwtp_name == "MWRDGC CALUMET WRP"), aes(x = sample_collect_date, y = voc_abundance, col = voc_abundance_NA)) + geom_point() 

#Interpolation 2: Electric Boogaloo
################################################################################
#combos2 <- mwrd_map_joined %>%
#  complete(wwtp_name, sample_collect_date)


grouped2 <- combos %>%
  group_by(wwtp_name, sample_collect_date)

grouped2$abvec <- replace_na_with(grouped2$abvec, 0)

no_samples <- grouped2 %>%
  summarize(voc_abundance = sum(abvec)) %>%
  mutate(wwtp_name = factor(wwtp_name)) %>%
  filter(voc_abundance <= 1 | is.na(voc_abundance)) %>% 
  mutate(sitedate = paste(wwtp_name, sample_collect_date)) %>%
  filter(voc_abundance != 0)

voc_sums2 <- voc_sums %>%
  mutate(sitedate = paste(wwtp_name, sample_collect_date))

voc_sums2[voc_sums2$sitedate %in% no_samples$sitedate & is.na(voc_sums2$voc_abundance),]$voc_abundance <- 0

filled2 <- voc_sums2 %>%
  bind_shadow() %>%
  group_by(named_variant_id, wwtp_name) %>%
  fill(voc_abundance, .direction = "down")

#Making the gifs
###############################################################################

filled_with_names <- filled2 %>%
  mutate(week = week(sample_collect_date)) %>%
  mutate(year = year(sample_collect_date)) %>%
  mutate(weeks_to_add = ((year - 2021) * 52)) %>%
  mutate(weeks_since_start = week + weeks_to_add)
filled_with_names$NAME <- recode(as.factor(filled2$wwtp_name), 
                                 "MWRDGC-KIRIE" = "Kirie",
                                 "MWRDGC HANOVER PARK WRP" = "Hanover",
                                 "MWRDGC STICKNEY WRP - SWS(N/S)" = "Stickney",
                                 "MWRDGC TERRENCE J O'BRIEN WTR RECLAMATION PLANT" = "OBrien",
                                 "MWRDGC-JOHN E. EGAN WRP" = "Egan",
                                 "MWRDGC LEMONT WRP" = "Lemont",
                                 "MWRDGC CALUMET WRP" = "Calumet",
                                 "MWRDGC STICKNEY WRP - WS" = "Stickney")

filled_joined <- right_join(filled_with_names, mwrd_map, by = "NAME") 

filled_joined$named_variant_id <- gsub("Omicron NA", "Omicron [other]", 
                                          filled_joined$named_variant_id)
filled_joined$named_variant_id <- gsub(" NA", "", filled_joined$named_variant_id) 
filled_joined$named_variant_id <- gsub("other", "Unnamed", filled_joined$named_variant_id) 
filled_joined$named_variant_id <- gsub("Other XBB", "Non-Omicron XBB", filled_joined$named_variant_id) 
filled_joined$named_variant_id <- gsub("Other recombinant", "Non-XBB recombinant", filled_joined$named_variant_id) 


chi_by_date <- filled_joined %>%
  ungroup() %>%
  arrange(sample_collect_date) %>%
  filter(voc_abundance > 0)

chi_order <- unique(chi_by_date$named_variant_id) %>% as.character()

filled_joined$named_variant_id <- factor(filled_joined$named_variant_id, 
                                            levels = order)
filled_joined <- filled_joined[!is.na(filled_joined$named_variant_id),]

all_vocs <- filled_joined %>%
  group_by(wwtp_name, named_variant_id) %>%
  mutate(rollav = rollmean(voc_abundance, k = 5, fill = 0, na.pad = TRUE, align = 'right')) %>%
  st_as_sf()

ggplot(filter(all_vocs, named_variant_id == "Omicron BA.5" & wwtp_name == "MWRDGC CALUMET WRP"), aes(x = sample_collect_date, y = voc_abundance)) + geom_point(aes(col = voc_abundance_NA))

dec142022 <- filter(all_vocs, sample_collect_date == "2022-12-14")

jan302022 <- filter(all_vocs, sample_collect_date == "2022-01-30")

chi_all_dates <- unique(all_vocs$sample_collect_date)
chi_dates01 <- all_dates[seq(1, length(chi_all_dates), 12)] 
chi_all_vocs_sample <- filter(all_vocs, sample_collect_date %in% chi_dates01)


ba5_anim = tm_shape(ba5_map) + tm_borders + tm_polygons(col = "rollav", palette = "viridis") +
  tm_facets(along = "weeks_since_start", free.coords = FALSE) #using "along =" instead of "by =" generates animation. "free.coords = FALSE" keeps the viewing window stationary.


tm_shape(jan302022) + tm_borders + tm_polygons(col = "rollav", palette = "viridis") +
  tm_facets(by = "named_variant_id")

all_map <- tm_shape(chi_all_vocs_sample) + tm_borders + tm_polygons(col = "rollav", palette = "viridis") +
  tm_facets(along = "sample_collect_date", by = "named_variant_id")

tmap_animation(all_map, filename = "mwrd_maps_sorted.gif", delay = 25)


ba5_map_grouped <- ba5_map %>%
  group_by(wwtp_name, weeks_since_start, sample_collect_date, named_variant_id)

ba5_anim = tm_shape(ba5_map) + tm_borders + tm_polygons(col = "rollav", palette = "viridis") +
  tm_facets(along = "sample_collect_date", free.coords = FALSE) #using "along =" instead of "by =" generates animation. "free.coords = FALSE" keeps the viewing window stationary.

ba5_anim

##############################################################

idph_ww_lin_data <- long_ww_lin_w_sample_info %>%
  right_join(idph_sf, by = "wwtp_name")

#Interpolation
###############################################################################
il_combos <- idph_ww_lin_data %>%
  complete(wwtp_name, sample_collect_date, named_variant_id)
il_grouped <- il_combos %>%
  group_by(wwtp_name, sample_collect_date, named_variant_id)

il_voc_sums <- il_grouped %>%
  summarize(voc_abundance = sum(abvec)) %>%
  mutate(named_variant_id = factor(named_variant_id)) %>%
  mutate(wwtp_name = factor(wwtp_name)) %>%
  filter(voc_abundance <= 1 | is.na(voc_abundance))

filled <- il_voc_sums %>%
  bind_shadow() %>%
  mutate(numdat = as.numeric(sample_collect_date)) %>%
  group_by(named_variant_id, wwtp_name) %>%
  fill(voc_abundance, .direction = "down")



###Interpolation 2
il_grouped2 <- il_combos %>%
  group_by(wwtp_name, sample_collect_date)

il_grouped2$abvec <- replace_na_with(il_grouped2$abvec, 0)

il_no_samples <- il_grouped2 %>%
  summarize(voc_abundance = sum(abvec)) %>%
  mutate(wwtp_name = factor(wwtp_name)) %>%
  filter(voc_abundance <= 1 | is.na(voc_abundance)) %>% 
  mutate(sitedate = paste(wwtp_name, sample_collect_date)) %>%
  filter(voc_abundance != 0)

il_voc_sums2 <- il_voc_sums %>%
  mutate(sitedate = paste(wwtp_name, sample_collect_date))

il_voc_sums2[il_voc_sums2$sitedate %in% il_no_samples$sitedate & is.na(il_voc_sums2$voc_abundance),]$voc_abundance <- 0

il_filled2 <- il_voc_sums2 %>%
  bind_shadow() %>%
  group_by(named_variant_id, wwtp_name) %>%
  fill(voc_abundance, .direction = "down")


ggplot(il_filled2, aes(x = sample_collect_date, y = voc_abundance)) + 
  geom_point(aes(col = named_variant_id, shape = voc_abundance_NA)) + facet_wrap(~wwtp_name)

####Maps
il_filled_with_names <- il_filled2 %>%
  mutate(week = week(sample_collect_date)) %>%
  mutate(year = year(sample_collect_date)) %>%
  mutate(weeks_to_add = ((year - 2021) * 52)) %>%
  mutate(weeks_since_start = week + weeks_to_add)


il_filled_joined <- left_join(il_filled_with_names, idph_sf, by = "wwtp_name") 



il_filled_joined$named_variant_id <- gsub("Omicron NA", "Omicron [other]", 
                                          il_filled_joined$named_variant_id)
il_filled_joined$named_variant_id <- gsub(" NA", "", il_filled_joined$named_variant_id) 
il_filled_joined$named_variant_id <- gsub("other", "Unnamed", il_filled_joined$named_variant_id) 
il_filled_joined$named_variant_id <- gsub("Other XBB", "Non-Omicron XBB", il_filled_joined$named_variant_id) 
il_filled_joined$named_variant_id <- gsub("Other recombinant", "Non-XBB recombinant", il_filled_joined$named_variant_id) 


by_date <- il_filled_joined %>%
  ungroup() %>%
  arrange(sample_collect_date) %>%
  filter(voc_abundance > 0)

order <- unique(by_date$named_variant_id) %>% as.character()

il_filled_joined$named_variant_id <- factor(il_filled_joined$named_variant_id, 
                                            levels = order)
il_filled_joined <- il_filled_joined[!is.na(il_filled_joined$named_variant_id),]

il_all_vocs <- il_filled_joined %>%
  group_by(wwtp_name, named_variant_id) %>%
  mutate(rollav = rollmean(voc_abundance, k = 5, fill = 0, na.pad = TRUE, align = 'right')) %>%
  st_as_sf()

ggplot(filter(il_all_vocs, named_variant_id == "Omicron BA.5" & wwtp_name == "MWRDGC CALUMET WRP"), aes(x = sample_collect_date, y = voc_abundance)) + geom_point(aes(col = voc_abundance_NA))

dec142022 <- filter(il_all_vocs, sample_collect_date == "2022-12-14")

jan302022 <- filter(il_all_vocs, sample_collect_date == "2022-01-30")

all_dates <- unique(il_all_vocs$sample_collect_date)
dates01 <- all_dates[seq(1, length(all_dates), 12)] 
il_all_vocs_sample <- filter(il_all_vocs, sample_collect_date %in% dates01)

tm_shape(IL) + tm_polygons(alpha = 0.2) +
  tm_shape(jan302022) + tm_borders + tm_dots(size = "rollav", col = "red") +
  tm_facets(by = "named_variant_id")

all_map <- tm_shape(all_vocs) + tm_borders + tm_polygons(col = "rollav", palette = "viridis") +
  tm_facets(along = "sample_collect_date", by = "named_variant_id")

il_all_map <- tm_shape(IL) + tm_polygons(alpha = 0.2) +
  tm_shape(il_all_vocs_sample) + tm_borders + tm_dots(size = "rollav", title = "Relative abundance", col = "red") +
  tm_facets(along = "sample_collect_date", by = "named_variant_id")

tmap_animation(il_all_map, filename = "bigger_sample_idph_anim_sorter.gif", delay = 25)