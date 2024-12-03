non_wrps <- c("MHOL23-00116", "SMH000072456", "SMH000111150", "SMH000116477", 
              "SMH000025183", "SMH000017048", "SMH000062398", "SMH000079342", 
              "CCJ - UNKNOWN ", "MHOL12-00024", "Schapers VA Home", 
              "Markword VA Home", "McAuley – Misericordia", "MHOL12-00024", 
              "MHOL23-00259", "CCJ- Division 5", "CCJ- Division 8 (RTU)", 
              "MHOL16-00342", "CCJ- Division 11 (Pod D)", "CCJ- Division 9", 
              "CCJ- Division 3 -Annex", "CCJ- Division 11 (Pod C)", 
              "CCJ- Division 2", "CCJ- Division 10", "CCJ- Division 6 (South)", 
              "", NA)


long_ww_lin_2 <- long_ww_lin_w_sample_info %>%
  filter(!is.na(sample_collect_date)) %>%
  mutate(dweek = as.numeric(sample_collect_date - min(sample_collect_date)) %/% 7)

ww_spread <- long_ww_lin_2 %>%
  filter(!(wwtp_name %in% non_wrps)) %>%
  group_by(aliases_removed, dweek) %>%
  summarize(number = length(unique(wwtp_name)))

voc_ww_spread <- long_ww_lin_2 %>%
  filter(!(wwtp_name %in% non_wrps)) %>%
  group_by(named_variant_id, dweek) %>%
  summarize(number = length(unique(wwtp_name)))

library(RColorBrewer)
colors = c(brewer.pal(n = 14, name="Set3"), brewer.pal(n = 8, name = "Set1"))
ggplot(voc_ww_spread, aes(x = dweek, y = rollmean(number, 3, na.pad = T), color = named_variant_id)) + geom_line(size = 1) +
  scale_color_manual(values = colors) +
  labs(color = "Variant name") +
  theme_minimal() +
  xlab("Weeks since Nov 8 2021") +
  ylab("Number of WWTPs detected")
