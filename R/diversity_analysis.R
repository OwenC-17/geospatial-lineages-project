library(tidyverse)
library(vegan)

non_wrps <- c("MHOL23-00116", "SMH000072456", "SMH000111150", "SMH000116477", 
              "SMH000025183", "SMH000017048", "SMH000062398", "SMH000079342", 
              "CCJ - UNKNOWN ", "MHOL12-00024", "Schapers VA Home", 
              "Markword VA Home", "McAuley – Misericordia", "MHOL12-00024", 
              "MHOL23-00259", "CCJ- Division 5", "CCJ- Division 8 (RTU)", 
              "MHOL16-00342", "CCJ- Division 11 (Pod D)", "CCJ- Division 9", 
              "CCJ- Division 3 -Annex", "CCJ- Division 11 (Pod C)", 
              "CCJ- Division 2", "CCJ- Division 10", "CCJ- Division 6 (South)", 
              "", NA)

ww_lin_data <- read_csv("../cleaned-formatted-data/sequencing-results-formatted-cleaned.csv", col_names = TRUE) 
#Later: must check sample_arrival_temp data not being erased. 

idph_ww_lin_data <- ww_lin_data %>%
  filter(sample_id != 19656) %>%
  right_join(idph_sf, by = "wwtp_name")

sample_diversity <- idph_ww_lin_data %>%
  group_by(sample_id) %>%
  summarise(richness = specnumber(abvec),
            shannon = diversity(abvec, index = "shannon"),
            simpson = diversity(abvec, index = "simpson"),
            inv_simpson = 1/simpson,
            Abundance = sum(abvec)) %>%
  left_join(idph_ww_lin_data, by = "sample_id")


diversity_summary <- sample_diversity %>%
  group_by(wwtp_name) %>%
  summarise(meanrichness = mean(richness),
            medrichness = median(richness),
            sdrichness = sd(richness),
            meanshannon = mean(shannon),
            medshannon = median(shannon),
            sdshannon = sd(shannon),
            pop = mean(population_served)) %>%
  arrange(desc(medrichness))

diversity_summary %>%
  left_join(leading_wwtps, by = "wwtp_name") %>%
  ggplot(aes(x = meanshannon, y = num_first)) + 
  geom_point() + 
  scale_y_log10() + 
  scale_x_log10()




sample_diversity$nice_wwtp_name <- factor(sample_diversity$nice_wwtp_name, 
                                          levels = richness_summary$nice_wwtp_name)

ggplot(sample_diversity, aes(x = nice_wwtp_name, y = richness)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggplot(sample_diversity, aes(x = nice_wwtp_name, y = shannon)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggplot(sample_diversity, aes(x = nice_wwtp_name, y = simpson)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  scale_y_log10()

ggplot(sample_diversity, aes(x = nice_wwtp_name, y = inv_simpson)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  scale_y_log10()

richness_summary %>%
  left_join(leading_wwtps, by = "wwtp_name") %>%
  ggplot(aes(x = meanrichness, y = num_first)) + 
  geom_point() + 
  scale_y_log10() + 
  scale_x_log10()
