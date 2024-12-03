first_chicago_ww_occurences_by_site <- long_ww_chicago %>%
  group_by(wwtp_name, aliases_removed) %>%
  filter(sample_collect_date == min(sample_collect_date))
