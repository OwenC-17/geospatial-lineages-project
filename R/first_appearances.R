library(tidyverse)
library(sf)
library(tmap)

ww_lin_data <- read_csv("../cleaned-formatted-data/sequencing-results-formatted-cleaned.csv", col_names = TRUE) 
#Later: must check sample_arrival_temp data not being erased. 

idph_ww_lin_data <- ww_lin_data %>%
  filter(sample_id != 19656) %>%
  right_join(idph_sf, by = "wwtp_name")

#Just a map of the locations without other data (for reference)
tm_shape(IL) + tm_polygons(alpha = 0.2) +
  tm_shape(leading_wwtps_no_singletons) + tm_dots(col = "cyan3", size = 0.5) +
  tm_layout(main.title = "WWTP locations") + tm_scale_bar(breaks = c(0, 50, 100), text.size = .8, position = c("left", "bottom"))


#NOTE: If a lineage appeared in two sites at once, both are counted!
first_appearances_by_lin <- idph_ww_lin_data %>%
  group_by(aliases_removed) %>%
  filter(weeks_since_start == min(weeks_since_start)) %>%
  distinct()

leading_wwtps <- first_appearances_by_lin %>%
  group_by(wwtp_name) %>%
  summarize(num_first = n()) %>%
  left_join(idph_sf, by = "wwtp_name") %>%
  st_as_sf()

tm_shape(IL) + tm_polygons(alpha = 0.2) +
  tm_shape(leading_wwtps) + tm_dots(col = "red", size = "num_first")

num_appearances <- ww_lin_data %>%
  group_by(aliases_removed) %>%
  summarize(nappearances = n())

singletons <- num_appearances %>%
  filter(nappearances == 1)

doubletons <- num_appearances %>%
  filter(nappearances <= 2)

tripletons <- num_appearances %>%
  filter(nappearances <= 3)

leading_wwtps_no_singletons <- first_appearances_by_lin %>%
  filter(!(aliases_removed %in% singletons$aliases_removed)) %>%
  group_by(wwtp_name) %>%
  summarize(num_first = n()) %>%
  left_join(idph_sf, by = "wwtp_name") %>%
  st_as_sf()

tm_shape(IL) + tm_polygons(alpha = 0.2) +
  tm_shape(leading_wwtps_no_singletons) + tm_dots(col = "red", size = "num_first")

######

leading_wwtps_no_doubletons <- first_appearances_by_lin %>%
  filter(!(aliases_removed %in% doubletons$aliases_removed)) %>%
  group_by(wwtp_name) %>%
  summarize(num_first = n()) %>%
  _join(idph_sf, by = "wwtp_name") %>%
  st_as_sf()

leading_wwtps_no_tripletons <- first_appearances_by_lin %>%
  filter(!(aliases_removed %in% tripletons$aliases_removed)) %>%
  group_by(wwtp_name) %>%
  summarize(num_first = n()) %>%
  left_join(idph_sf, by = "wwtp_name") %>%
  st_as_sf()

##plots
################

all_lineages = tm_shape(IL) + tm_polygons(alpha = 0.2) +
  tm_shape(leading_wwtps) + tm_dots(col = "population_served", palette = "viridis", size = "num_first") +
  tm_layout(main.title = "No restriction")

no_singletons = tm_shape(IL) + tm_polygons(alpha = 0.2) +
  tm_shape(leading_wwtps_no_singletons) + tm_dots(col = "population_served", palette = "viridis", size = "num_first") +
  tm_layout(main.title = "> 1 appearances")

no_doubletons = tm_shape(IL) + tm_polygons(alpha = 0.2) +
  tm_shape(leading_wwtps_no_doubletons) + tm_dots(col = "population_served", palette = "viridis", size = "num_first") +
  tm_layout(main.title = "> 2 appearances")

no_tripletons = tm_shape(IL) + tm_polygons(alpha = 0.2) +
  tm_shape(leading_wwtps_no_tripletons) + tm_dots(col = "population_served", palette = "viridis", size = "num_first") +
  tm_layout(main.title = "> 3 appearances")

tmap_arrange(all_lineages, no_singletons, no_doubletons)

ggplot(leading_wwtps, aes(x = population_served, y = num_first)) + geom_point() + geom_smooth(method = "lm") +
  scale_x_log10() +
  scale_y_log10() +
  xlab("population served") +
  ylab("number of first detections") +
  theme_bw()

pop_first_model <- lm(log(num_first) ~ log(population_served), data = leading_wwtps)
summary(pop_first_model)
plot(pop_first_model)

################################################################################
###Looking at Named Variants (VOCs, VOIs, & VUMs)
first_appearances_by_voc <- idph_ww_lin_data %>%
  select(named_variant_id, weeks_since_start, wwtp_name, nice_wwtp_name, capacity_mgd) %>%
  distinct() %>%
  group_by(named_variant_id) %>%
  filter(weeks_since_start == min(weeks_since_start)) %>%
  left_join(idph_sf, by = "wwtp_name") %>%
  st_as_sf()

leading_wwtps_voc <- first_appearances_by_voc %>%
  group_by(wwtp_name) %>%
  summarize(num_first = n()) %>%
  right_join(as.data.frame(idph_sf, by = "wwtp_name")) %>%
  st_as_sf()

tm_shape(IL) + tm_polygons(alpha = 0.2) +
  tm_shape(leading_wwtps_voc) + tm_dots(col = "log10_pop", breaks = c(3.0, 4.0, 5.0, 6.0, 7.0), palette = "viridis", size = "num_first") +
  tm_layout(main.title = "Grouped to \nNamed Variants")

################################################################################
#Repeat analyses with different tax levels:


#Top tier#########################################
first_appearances_by_tier1 <- idph_ww_lin_data %>%
  select(top_lin_id, weeks_since_start, wwtp_name, nice_wwtp_name, capacity_mgd) %>%
  distinct() %>%
  group_by(top_lin_id) %>%
  filter(weeks_since_start == min(weeks_since_start)) %>%
  left_join(idph_sf, by = "wwtp_name") %>%
  st_as_sf()

leading_wwtps_tier1 <- first_appearances_by_tier1 %>%
  group_by(wwtp_name) %>%
  summarize(num_first = n()) %>%
  right_join(as.data.frame(idph_sf, by = "wwtp_name")) %>%
  st_as_sf()

tm_shape(IL) + tm_polygons(alpha = 0.2) +
  tm_shape(leading_wwtps_tier1) + tm_dots(col = "log10_pop", breaks = c(3.0, 4.0, 5.0, 6.0, 7.0), palette = "viridis", size = "num_first") +
  tm_layout(main.title = "Grouped to Tier 1")


#Second tier############################################

idph_ww_lin_data <- idph_ww_lin_data %>%
  unite("top_lin_id", "sub1", col = "tier2id", sep = ".", remove = FALSE, na.rm = TRUE)

first_appearances_by_tier2 <- idph_ww_lin_data %>%
  select(tier2id, weeks_since_start, wwtp_name, nice_wwtp_name, capacity_mgd) %>%
  distinct() %>%
  group_by(tier2id) %>%
  filter(weeks_since_start == min(weeks_since_start)) %>%
  left_join(idph_sf, by = "wwtp_name") %>%
  st_as_sf()

leading_wwtps_tier2 <- first_appearances_by_tier2 %>%
  group_by(wwtp_name) %>%
  summarize(num_first = n()) %>%
  right_join(as.data.frame(idph_sf, by = "wwtp_name")) %>%
  st_as_sf()

tm_shape(IL) + tm_polygons(alpha = 0.2) +
  tm_shape(leading_wwtps_tier2) + tm_dots(col = "log10_pop", breaks = c(3.0, 4.0, 5.0, 6.0, 7.0), palette = "viridis", size = "num_first") +
  tm_layout(main.title = "Grouped to Tier 2")


################################################################################
#Examine WRP-wise sampling intensity:
first_loc_date <- ww_lin_data %>%
  group_by(wwtp_name) %>%
  summarize(first_date = min(sample_collect_date))


