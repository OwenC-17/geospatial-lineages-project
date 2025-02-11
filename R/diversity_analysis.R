library(tidyverse)
library(sf)
library(ade4)
library(vegan)
library(gclus)
library(cluster)
library(FD)
library(RColorBrewer)
library(paletteer)
library(zoo)

non_wrps <- c("MHOL23-00116", "SMH000072456", "SMH000111150", "SMH000116477",
              "SMH000025183", "SMH000017048", "SMH000062398", "SMH000079342", 
              "CCJ - UNKNOWN ", "MHOL12-00024", "Schapers VA Home", 
              "Markword VA Home", "McAuley – Misericordia", "MHOL12-00024", 
              "MHOL23-00259", "CCJ- Division 5", "CCJ- Division 8 (RTU)", 
              "MHOL16-00342", "CCJ- Division 11 (Pod D)", "CCJ- Division 9", 
              "CCJ- Division 3 -Annex", "CCJ- Division 11 (Pod C)", 
              "CCJ- Division 2", "CCJ- Division 10", "CCJ- Division 6 (South)", 
              "", NA)

ww_lin_data <- read_csv("../cleaned-formatted-data/sequencing-results-formatted-cleaned.csv", col_names = TRUE, guess_max = 1000000) 
#Later: must check sample_arrival_temp data not being erased. 
ww_lin_data <- enforce_types(ww_lin_data)

wwtp_names_dictionary <- read_csv("../input/wwtp_names_dictionary.csv")

ww_lin_data$site_type <- wwtp_names_dictionary$site_type[match(
  ww_lin_data$wwtp_name, wwtp_names_dictionary$wwtp_name)] %>%
  factor(ordered = FALSE)

idph_ww_lin_data <- ww_lin_data %>%
  filter(sample_id != 19656 & !is.na(sample_id)) %>%
  right_join(idph_sf, by = "wwtp_name") %>%
  distinct()

sample_diversity <- idph_ww_lin_data %>%
  group_by(sample_id) %>%
  summarise(richness = specnumber(abvec),
            shannon = diversity(abvec, index = "shannon"),
            simpson = diversity(abvec, index = "simpson"),
            inv_simpson = 1/simpson,
            Abundance = sum(abvec)) %>%
  left_join(idph_ww_lin_data, by = "sample_id")


diversity_summary <- sample_diversity %>%
  group_by(nice_wwtp_name) %>%
  summarise(meanrichness = mean(richness),
            medrichness = median(richness),
            sdrichness = sd(richness),
            meanshannon = mean(shannon),
            medshannon = median(shannon),
            sdshannon = sd(shannon),
            pop = mean(population_served)) %>%
  arrange(desc(medrichness))

sample_diversity$nice_wwtp_name <- factor(sample_diversity$nice_wwtp_name, levels = diversity_summary$nice_wwtp_name)

ggplot(sample_diversity, aes(x = nice_wwtp_name, y = richness)) + 
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_y_log10()

ggplot(sample_diversity, aes(x = nice_wwtp_name, y = shannon)) + 
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_y_log10()

ggplot(sample_diversity, aes(x = nice_wwtp_name, y = simpson)) + 
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_y_log10()

ggplot(sample_diversity, aes(x = nice_wwtp_name, y = inv_simpson)) + 
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_y_log10()

ggplot(sample_diversity, aes(x = sample_type, y = richness)) + 
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

ggplot(sample_diversity, aes(x = sample_type, y = shannon)) + 
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

ggplot(sample_diversity, aes(x = sample_collect_date, y = shannon)) + geom_point() + geom_smooth()

ggplot(sample_diversity, aes(x = sample_collect_date, y = shannon)) + 
  geom_point(aes(col = log(population_served))) +
  scale_color_paletteer_d(palette = "viridis::viridis")

ggplot(sample_diversity, aes(x = sample_collect_date, y = shannon)) + 
  geom_line(aes(col = nice_wwtp_name, alpha = 0.5)) +
  theme(legend.position = "none")
ggplot(sample_diversity, aes(x = sample_collect_date, y = rollmean(shannon, 30, na.pad = TRUE))) + 
  geom_line(aes(col = nice_wwtp_name, alpha = 0.5)) +
  theme(legend.position = "none")


ggplot(sample_diversity, aes(x = sample_collect_date, y = log(richness))) + 
  geom_point(aes(col = log(population_served))) +
  scale_color_paletteer_c(palette = "viridis::viridis")


ggplot(sample_diversity, aes(x = sample_collect_date, y = richness)) + geom_smooth()



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

################################################################################

idph_by_week <- idph_ww_lin_data %>%
  split(idph_ww_lin_data$weeks_since_start)

SAMPLEWEEK <- idph_by_week[[73]]

WIDESAMPLE <- SAMPLEWEEK %>%
  select(-c(linvec, full_lineage_id, top_lin_id, sub1, sub2, sub3, named_variant_id)) %>%
  pivot_wider(names_from = aliases_removed, values_from = abvec)

ABUNDANCES_ONLY <- WIDESAMPLE %>%
  select(sample_id, XBB.1.9.2.1:XBB.1.9.1.4)
ABUNDANCES_ONLY[is.na(ABUNDANCES_ONLY)] <- 0
rownames(ABUNDANCES_ONLY) <- ABUNDANCES_ONLY$sample_id
ABUND_MATRIX <- as.matrix(ABUNDANCES_ONLY[,2:ncol(ABUNDANCES_ONLY)])
rownames(ABUND_MATRIX) <- ABUNDANCES_ONLY$sample_id

BRAY <- vegdist(ABUND_MATRIX)
LOGBRAY <- vegdist(log1p(ABUND_MATRIX))

PCOA_BRAY <- cmdscale(BRAY) %>% as.data.frame()
PCOA_BRAY$sample_id <- rownames(PCOA_BRAY)

JOINED_PCOA <- PCOA_BRAY %>%
  left_join(WIDESAMPLE, by = "sample_id")

ggplot(JOINED_PCOA, aes(x = V1, y = V2)) + geom_point(aes(colour = sample_type))

PCOA_LOGBRAY <- cmdscale(LOGBRAY) %>% as.data.frame()
PCOA_LOGBRAY$sample_id <- rownames(PCOA_LOGBRAY)

JOINED_PCOALOG <- PCOA_LOGBRAY %>%
  left_join(WIDESAMPLE, by = "sample_id")

ggplot(JOINED_PCOALOG, aes(x = V1, y = V2)) + geom_point(aes(colour = sample_type))

SCALED <- scale(ABUNDANCES_ONLY[,2:916])

#All of the above on all dates
################################################################################

WIDEALL <- idph_ww_lin_data %>%
  select(-c(linvec, full_lineage_id, top_lin_id, sub1, sub2, sub3, named_variant_id)) %>%
  pivot_wider(names_from = aliases_removed, values_from = abvec, values_fn = mean)

ALL_ABUNDANCES_ONLY <- WIDEALL %>%
  select(sample_id, B.1.1.529.1.1.9:XDZ) %>%
  filter(!is.na(sample_id))
ALL_ABUNDANCES_ONLY[is.na(ALL_ABUNDANCES_ONLY)] <- 0
rownames(ALL_ABUNDANCES_ONLY) <- ALL_ABUNDANCES_ONLY$sample_id
ALL_ABUND_MATRIX <- as.matrix(ALL_ABUNDANCES_ONLY[,2:ncol(ALL_ABUNDANCES_ONLY)])
rownames(ALL_ABUND_MATRIX) <- ALL_ABUNDANCES_ONLY$sample_id

BIGBRAY <- vegdist(ALL_ABUND_MATRIX)
BIGLOGBRAY <- vegdist(log1p(ALL_ABUND_MATRIX))

PCOA_BIGBRAY <- cmdscale(BIGBRAY) %>% as.data.frame()
PCOA_BIGBRAY$sample_id <- rownames(PCOA_BIGBRAY)

JOINED_BIGPCOA <- PCOA_BIGBRAY %>%
  left_join(WIDEALL, by = "sample_id")

ggplot(JOINED_BIGPCOA, aes(x = V1, y = V2)) + geom_point(aes(colour = nice_wwtp_name))

PCOA_BIGLOGBRAY <- cmdscale(BIGLOGBRAY) %>% as.data.frame()
PCOA_BIGLOGBRAY$sample_id <- rownames(PCOA_BIGLOGBRAY)

JOINED_BIGPCOALOG <- PCOA_BIGLOGBRAY %>%
  left_join(WIDEALL, by = "sample_id")

ggplot(JOINED_BIGPCOALOG, aes(x = V1, y = V2)) + geom_point(aes(colour = sample_type))

SCALED <- scale(ABUNDANCES_ONLY[,2:916])
