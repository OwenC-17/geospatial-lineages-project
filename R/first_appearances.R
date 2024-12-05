ww_lin_data <- read_csv("../../ClinicalVsWW/ww_lineages_and_metadata_complete.csv",col_names = TRUE) %>%
  mutate_at(.vars = "sample_arrival_temp", .funs = gsub, #mutate_at allows selecting a column instead of making a new one 
            pattern = "\\.{2}|,", replacement = "\\.")  %>%
  mutate(weeks_to_add = ((year - 2021) * 52)) %>%
  mutate(weeks_since_start = week + weeks_to_add)

ww_lin_data2 <- read_csv("../cleaned-formatted-data/sequencing-results-formatted-cleaned.csv",col_names = TRUE) %>%
  mutate_at(.vars = "sample_arrival_temp", .funs = gsub, #mutate_at allows selecting a column instead of making a new one 
            pattern = "\\.{2}|,", replacement = "\\.")  %>%
  mutate(weeks_to_add = ((year - 2021) * 52)) %>%
  mutate(weeks_since_start = week + weeks_to_add)


idph_ww_lin_data <- ww_lin_data %>%
  filter(sample_id != 19656) %>%
  right_join(idph_sf, by = "wwtp_name")

first_appearances_by_lin <- idph_ww_lin_data %>%
  group_by(long_lineage_id) %>%
  filter(sample_collect_date == min(sample_collect_date)) %>%
  distinct()

leading_wwtps <- first_appearances_by_lin %>%
  group_by(wwtp_name) %>%
  summarize(num_first = n()) %>%
  left_join(idph_sf, by = "wwtp_name") %>%
  st_as_sf()

tm_shape(IL) + tm_polygons(alpha = 0.2) +
  tm_shape(leading_wwtps) + tm_dots(col = "red", size = "num_first")

num_appearances <- first_appearances_by_lin %>%
  group_by(long_lineage_id) %>%
  summarize(nappearances = n())

singletons <- num_appearances %>%
  filter(nappearances == 1)

doubletons <- num_appearances %>%
  filter(nappearances <= 2)

tripletons <- num_appearances %>%
  filter(nappearances <= 3)

leading_wwtps_no_singletons <- first_appearances_by_lin %>%
  filter(!(long_lineage_id %in% singletons$long_lineage_id)) %>%
  group_by(wwtp_name) %>%
  summarize(num_first = n()) %>%
  left_join(idph_sf, by = "wwtp_name") %>%
  st_as_sf()

tm_shape(IL) + tm_polygons(alpha = 0.2) +
  tm_shape(leading_wwtps_no_singletons) + tm_dots(col = "red", size = "num_first")


tm_shape(IL) + tm_polygons(alpha = 0.2) +
  tm_shape(leading_wwtps_no_singletons) + tm_dots(col = "cyan3", size = 0.5) +
  tm_layout(main.title = "WWTP locations") + tm_scale_bar(breaks = c(0, 50, 100), text.size = .8, position = c("left", "bottom"))

######

leading_wwtps_no_doubletons <- first_appearances_by_lin %>%
  filter(!(long_lineage_id %in% doubletons$long_lineage_id)) %>%
  group_by(wwtp_name) %>%
  summarize(num_first = n()) %>%
  left_join(idph_sf, by = "wwtp_name") %>%
  st_as_sf()

leading_wwtps_no_tripletons <- first_appearances_by_lin %>%
  filter(!(long_lineage_id %in% tripletons$long_lineage_id)) %>%
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

first_appearances_by_voc <- idph_ww_lin_data %>%
  group_by(voc_id) %>%
  filter(sample_collect_date == min(sample_collect_date)) %>%
  group_by(voc_id, wwtp_name, sample_collect_date) %>%
  summarize(voc_abund = sum(abvec)) %>%
  left_join(idph_sf, by = "wwtp_name") %>%
  st_as_sf()

leading_wwtps_voc <- first_appearances_by_voc %>%
  group_by(wwtp_name) %>%
  summarize(num_first = n()) %>%
  left_join(as.data.frame(idph_sf, by = "wwtp_name")) %>%
  st_as_sf()

tm_shape(IL) + tm_polygons(alpha = 0.2) +
  tm_shape(leading_wwtps_voc) + tm_dots(col = "population_served", palette = "viridis", size = "num_first") +
  tm_layout(main.title = "Grouped to VOCs")


tm_shape(IL) + tm_polygons(alpha = 0.2) +
  tm_shape(first_appearances_by_voc) + tm_dots(col = "population_served", palette = "viridis", size = 1) +
  tm_facets(along = "voc_id", free.coords = FALSE) +
  tm_layout(main.title = "Locations of first appearance") +
  tm_text(text = "sample_collect_date")


################################################################################
first_loc_date <- ww_lin_data %>%
  group_by(wwtp_name) %>%
  summarize(first_date = min(sample_collect_date))

