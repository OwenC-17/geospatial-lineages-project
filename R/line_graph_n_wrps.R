###This script produces a line graph of lineage prevalence over time.
###Run load_data.R in the current environment first. 
###Prevalence is here defined as the number of WRPs at which a 
###lineage was detected.

library(tidyverse)
library(zoo)

#We only want to visualize WRPs, so remove these.
#(This list was last updated for May 2024 data).
non_wrps <- c("MHOL23-00116", "SMH000072456", "SMH000111150", "SMH000116477", 
              "SMH000025183", "SMH000017048", "SMH000062398", "SMH000079342", 
              "CCJ - UNKNOWN ", "MHOL12-00024", "Schapers VA Home", 
              "Markword VA Home", "McAuley – Misericordia", "MHOL12-00024", 
              "MHOL23-00259", "CCJ- Division 5", "CCJ- Division 8 (RTU)", 
              "MHOL16-00342", "CCJ- Division 11 (Pod D)", "CCJ- Division 9", 
              "CCJ- Division 3 -Annex", "CCJ- Division 11 (Pod C)", 
              "CCJ- Division 2", "CCJ- Division 10", "CCJ- Division 6 (South)", 
              "", NA)

wwtp_names_dictionary <- read_csv("../input/wwtp_names_dictionary.csv")

long_ww_lin_w_sample_info$site_type <- wwtp_names_dictionary$site_type[match(
  long_ww_lin_w_sample_info$wwtp_name, wwtp_names_dictionary$wwtp_name)] %>%
  factor(ordered = FALSE)

##Note: dweek = number of weeks since start of sampling (2021-11-08), while
##      weeks_since_start is the number of weeks since 2021-01-01.

long_ww_lin_2 <- long_ww_lin_w_sample_info %>%
  filter(!is.na(sample_collect_date))

ww_spread <- long_ww_lin_2 %>%
  filter(site_type == "wwtp") %>%
  group_by(aliases_removed, weeks_since_start) %>%
  summarize(number = length(unique(wwtp_name)))

voc_ww_spread <- long_ww_lin_2 %>%
  filter(site_type == "wwtp") %>%
  group_by(named_variant_id, weeks_since_start) %>%
  summarize(number = length(unique(wwtp_name)))

library(RColorBrewer)
colors = c(brewer.pal(n = 12, name="Set3"), brewer.pal(n = 8, name = "Set2"), brewer.pal(n = 9, name = "Set1"), brewer.pal(n = 12, name = "Set3"), brewer.pal(n = 8, name = "Set2"))
ggplot(ww_spread, aes(x = weeks_since_start, y = rollmean(number, 3, na.pad = T), color = aliases_removed)) + geom_line(size = 1) +
  labs(color = "Variant name") +
  theme_minimal() +
  xlab("Weeks since Nov 8 2021") +
  ylab("Number of WWTPs detected")
