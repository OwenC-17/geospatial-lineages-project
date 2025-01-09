###NOTE: This script is for loading data that has already been formatted with
###import-data.R and saved as .csv. Do not add any other functions to this file. 

library(tidyverse)
library(lubridate)

long_ww_lin_w_sample_info <- read_csv(
  "../cleaned-formatted-data/sequencing-results-formatted-cleaned.csv",
  guess_max = 1000000) #Must be longer than the number of rows or some data get
                       #deleted due to incorrect type assignment

enforce_types <- function(long_ww_lin_w_sample_info) {
  long_ww_lin_w_sample_info <- long_ww_lin_w_sample_info %>%
    
    mutate(across(c("abvec", "coverage", "flow_rate",
                    "sample_arrival_temp",), as.double),
           
           across(c("linvec", "sample_id", "site_id", "wwtp_name", 
                    "sample_type", "notes", "nice_wwtp_name", "full_lineage_id",
                    "top_lin_id", "sub1", "sub2", "sub3", "aliases_removed",
                    "named_variant_id"), as.character),
           
           across(c("sample_collect_time"), hms),
           
           across(c("sample_collect_date", "sample_processed_date",
                    "test_result_date", "sample_received_date"), as_date),
           
           across(c("sample_collect_datetime"), as_datetime, tz = "US/Central"),
           
           across(c("year", "week", "days_since_start",
                    "weeks_since_start"), as.integer)
           )
  return(long_ww_lin_w_sample_info)
}

long_ww_lin_w_sample_info <- enforce_types(long_ww_lin_w_sample_info)
