library(tidyverse)
library(lubridate)

long_ww_lin_w_sample_info <- read.csv(
  "../cleaned-formatted-data/sequencing-results-formatted-cleaned.csv")

enforce_types <- function(long_ww_lin_w_sample_info) {
  long_ww_lin_w_sample_info <- long_ww_lin_w_sample_info %>%
    mutate(across(c("abvec", "coverage", "flow_rate",
                    "sample_arrival_temp",), as.double),
           across(c("wwtp_name", "nice_wwtp_name", "aliases_removed",
                     "sample_id", "site_id", "sample_type", "notes",
                     "full_lineage_id", "named_variant_id"), as.character),
           across(c("sample_collect_time"), hms),
           across(c("test_result_date", "sample_collect_date",
                     "sample_processed_date", "sample_received_date"), as_date))
  return(long_ww_lin_w_sample_info)
}
long_ww_lin_w_sample_info <- enforce_types(long_ww_lin_w_sample_info)

lineage_components <- (str_split(long_ww_lin_w_sample_info$full_lineage_id, "\\.",simplify = TRUE))

lineage_components[lineage_components == ""] <- NA

lineage_components <- base::apply(data.frame(lineage_components), 2, FUN = str_replace, pattern = "^$", replacement = NA)

long_ww_lin_w_sample_info <- long_ww_lin_w_sample_info %>%
  mutate(top_lin_id = lineage_components[,1],
         sub1 = lineage_components[,2],
         sub2 = lineage_components[,3],
         sub3 = lineage_components[,4])


#
mini <- head(long_ww_lin_w_sample_info[,24:27], 100)

#This is just to confirm that the separate lineages combine properly (should
#return 0).
sum(unite(long_ww_lin_w_sample_info[,24:27], 
          "combined", sep = ".", na.rm = TRUE)
    != long_ww_lin_w_sample_info$full_lineage_id)

