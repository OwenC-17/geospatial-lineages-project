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


