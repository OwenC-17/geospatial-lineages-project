chicago_sites <- read.csv("../input/chicago_sites.csv", header = FALSE, col.names = c("number", "sites")) %>%
  mutate(sites = gsub('\\"', '', sites))

long_ww_chicago <- filter(long_ww_lin_w_sample_info, 
                          wwtp_name %in% chicago_sites$sites) %>%
  mutate(wwtp_name = gsub('METROPOLITAN ', 'M', wwtp_name))

#list_of_others <- as.data.frame(table(filter(long_ww_lin_w_sample_info, named_variant_id == "other")$full_lineage_id))

