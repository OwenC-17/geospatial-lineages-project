###This script is used to import raw data (downloaded from the LIMS), clean, and
###format it. Any modifications to the data that take place BEFORE analysis go
###here. 

library("tidyverse")
library("zoo")
library("lubridate")

#Read in the data
################################################################################
unf_ww <- read_csv("../input/2025-01-03-lineage-results.csv",
                   col_types = "fffDdcccc--")


#Function to convert lineage names and abundances from lists to long
#tabular format:
################################################################################
tablify_lineage_abundances <- function(dfrow) {
  require("tidyverse")
  linvec <- gsub("\\[|]|\\s|\\'", "", dfrow$lineages) %>%
    gsub(",like[,]?(?=[0-9])", "~", ., perl = TRUE) %>%
    gsub(",like", "~", .) %>%
    strsplit(",") %>% #Split by commas, turn the string into a list
    unlist() #Turn the list into a character vector

  abvec <- gsub("\\[|]|\\s|\\'", "", dfrow$abundances) %>% #Same as above
    strsplit(",") %>%
    unlist() %>%
    as.numeric() #It's text by default, we need floats
  if (length(linvec) == length(abvec)) {
    df <- data.frame(linvec, abvec) %>% #Join lineages and abundances
      #Get IDs and coverage from original input:
      mutate(sample_id = dfrow$sample_id, coverage = dfrow$coverage)
    return(df)
  } else {
    message(row, ": ", linvec, " & ", abvec, " do not match")
  }
}

#Apply the function:
if (exists("long_ww_lin")) {
  #Have to delete this or else the for loop will just make it longer:
  rm(long_ww_lin)
}

for (row in seq_len(nrow(unf_ww))) { #First row is 1, not 0
  if (!exists("long_ww_lin")) { #Create df for output if needed
    long_ww_lin <<- tablify_lineage_abundances(unf_ww[row, ])
  } else {
    temporary <- tablify_lineage_abundances(unf_ww[row, ])
    long_ww_lin <<- bind_rows(long_ww_lin, temporary)
    #Join lengthened row to existing data
    rm(temporary) #Save memory <3
  }
}

rm(unf_ww) #This is big and we don't need it anymore

###Add metadata
################################################################################

#Read metadata file
#(This file associates sample IDs to all the data besides lineage abundances,
#including collection site, date, sample type, flow rate, temp etc.):

#Old metadata:
sample_metadata <- read_csv("../input/sample_metadata.csv",
                            col_types = "ff--ff-dDtD--Dccc--")

#New metadata (will generate a warning about an improperly formatted date, which
#is resolved next):
sample_metadata2 <- read_csv("../input/sample_metadata_20250103.csv",
                             col_types = "ff--ff-dDtD--Dccc--") %>%
  filter(!is.na(sample_id))

#These dates had typos that need to be fixed:
sample_metadata2[(sample_metadata2$sample_id == "102619"), ]$sample_collect_date <- as.Date("2024-08-20")
sample_metadata2[(sample_metadata2$sample_id == "105230"), ]$sample_collect_date <- as.Date("2024-10-16")

#Format raw metadata:
sample_metadata <- sample_metadata %>%
  mutate(sample_received_date = as.Date(sample_received_date)) %>%
  mutate(sample_arrival_temp = gsub("\\.\\.|,", "\\.", sample_arrival_temp)) %>%
  mutate(sample_arrival_temp = gsub("`|\\.$", "", sample_arrival_temp)) %>%
  mutate(sample_arrival_temp = as.double(sample_arrival_temp)) %>%
  mutate(year = year(sample_collect_date)) %>%
  mutate(week = week(sample_collect_date)) %>%
  mutate(weeks_to_add = ((year - 2021) * 52)) %>%
  mutate(weeks_since_start = week + weeks_to_add)


#Some values in sample_arrival_temp have issues preventing conversion to
#numeric values. Find and fix the issues:

#Get unique temperature values:
all_temps <- unique(sample_metadata2$sample_arrival_temp)
#Find values that don't properly convert to numeric:
weird_temps <- all_temps[is.na(as.numeric(all_temps))]
weird_temps

#Any string of consecutive dots or commas replaced with a single dot:
weird2 <- str_replace_all(weird_temps, "\\.+|,+", ".")
#For testing purposes, to be sure negative numbers won't be excluded:
weird2[2] <- "-12.3." 
weird2

#Any characters that are not numeric or dots or minus signs removed
#AND trailing dots removed:
weird3 <- str_remove_all(weird2, "[^[:digit:]\\.-]|(\\.$)")
weird3

#Compare original weird temps to corrected values to make sure they behave as 
#expected:
data.frame(weird_temps, as.numeric(weird3))

#Now fix the temperatures in the main data and make them numeric
#(Note use of nested pipes!):
sample_metadata2 <- sample_metadata2 %>%
  mutate(sample_arrival_temp = sample_arrival_temp %>%
           str_replace_all(pattern = "\\.+|,+",
                           replacement = ".") %>%
           str_remove_all("[^[:digit:]\\.-]|(\\.$)") %>%
           as.double())

#Create helpful alternative representations of dates:
sample_metadata2 <- sample_metadata2 %>%
  mutate(sample_received_date = as.Date(sample_received_date)) %>%
  mutate(year = year(sample_collect_date)) %>%
  filter(year != 1999) %>% #This will also remove samples with NA dates!
  mutate(week = week(sample_collect_date)) %>%
  mutate(days_since_start = sample_collect_date - min(sample_collect_date)) %>%
  mutate(weeks_since_start = difftime(sample_collect_date,
                                      min(sample_collect_date), 
                                      units = "weeks") %>% floor()) %>%
  #In the following line, a warning will be caused by samples without collect
  #times. This is fine.
  mutate(sample_collect_datetime = paste(sample_collect_date, 
                                         sample_collect_time) %>%
           ymd_hms())

#Merge redundant sample types (these exist due to inconsistent manual entries
#into the LIMS):
sample_metadata2$sample_type <- case_match(
  sample_metadata2$sample_type,
  c("Grab sample", "grab", "Grab") ~ "grab",
  c("24-hr flow-weighted composite", "Flow-proportional 24-hr composite",
    "Volume-proportional 24-hr composite") ~ "flow_weighted_composite",
  c("Time-proportional 24-hr composite", "Time proportional 24 Hr composite",
    "24-hr time-weighted composite") ~ "time_weighted_composite",
  "24-hr manual composite" ~ "manual_composite",
  c("Moore swab", "Moore Swab") ~ "moore_swab"
) %>% factor(ordered = FALSE)

#Merge redundant location names and make them nicer to read
wwtp_names_dictionary <- read_csv("../input/wwtp_names_dictionary.csv")
sample_metadata2$nice_wwtp_name <- wwtp_names_dictionary$pretty_name[match(
  sample_metadata2$wwtp_name, wwtp_names_dictionary$wwtp_name)] %>%
  factor(ordered = FALSE)

#Add the metadata to the lineage abundance table generated in previous
#section:
long_ww_lin_w_sample_info <- left_join(long_ww_lin,
                                       sample_metadata2,
                                       by = "sample_id")

#Cleanup
rm(long_ww_lin, sample_metadata2, wwtp_names_dictionary)


###Format lineage names
################################################################################

#Make a column specifically for unedited lineages (they will be edited soon)
long_ww_lin_w_sample_info$full_lineage_id <- long_ww_lin_w_sample_info$linvec


#Not currently using (use str_split to get lineage_components instead):
#Separate lineages by dots (this will generate a warning. It is okay.)
#long_ww_lin_w_sample_info <- long_ww_lin_w_sample_info %>%
#  separate(full_lineage_id,
#           paste0("name",
#                  seq_len(str_count(long_ww_lin_w_sample_info$full_lineage_id, "\\.")
#                          + 1)),
#           sep = "\\.", remove = FALSE)

lineage_components <- (str_split(long_ww_lin_w_sample_info$full_lineage_id, "\\.",simplify = TRUE))

lineage_components[lineage_components == ""] <- NA

lineage_components <- base::apply(data.frame(lineage_components), 2, FUN = str_replace, pattern = "^$", replacement = as.character(NA))

long_ww_lin_w_sample_info <- long_ww_lin_w_sample_info %>%
  mutate(top_lin_id = lineage_components[,1],
         sub1 = lineage_components[,2],
         sub2 = lineage_components[,3],
         sub3 = lineage_components[,4])

#This is just to confirm that the separate lineages combine properly (should
#return 0).
sum(unite(long_ww_lin_w_sample_info[,23:26], 
          "combined", sep = ".", na.rm = TRUE)
    != long_ww_lin_w_sample_info$full_lineage_id)


#Not using but may need for other functions, keep until confirmed unneeded:
#Create column with lineages only up to fist nested number:
#long_ww_lin_w_sample_info$t2lin <- paste(long_ww_lin_w_sample_info$name1,
#                                         long_ww_lin_w_sample_info$name2,
#                                         "x",
#                                         sep = ".")

###Some lineages have more than one name. The following code makes them
###consistent.

#Lineage notes from cov-lineages/pango-designation on Github (ignore parsing warning):
lineage_notes <- read_table("../input/lineage_notes.txt", skip = 1,
                            col_names = FALSE)

#Pull designated aliases from the lineage notes:
aliases <- lineage_notes %>%
  filter(X2 == "Alias" & X3 == "of") %>%
  select(all_of(c("X1", "X4"))) %>%
  rename(alias = "X1", original = "X4") %>%
  mutate(original = gsub(",", "", original))
rm(lineage_notes)

#Convert alias list to named character vector (faster than the DF):
alias_vector <- as.character(aliases$original)
names(alias_vector) <- aliases$alias

#Function to replace aliases in lineage IDs:
replace_alias <- function(id, a = alias_vector) {
  require("tidyverse")
  if (id %in% names(a)) {
    new_id <- a[id]
  } else {
    new_id <- id
  }
  message(id, "->", new_id)
  return(new_id)
}

#Use the above function to get consistent lineage IDs in the main data (a
#new data.frame is created for safety):
long_ww_lin_3 <- long_ww_lin_w_sample_info %>%
  mutate(aliases_removed = lapply(.$full_lineage_id, replace_alias))
long_ww_lin_3$aliases_removed <- unlist(long_ww_lin_3$aliases_removed)

#If above chunk runs successfully, modify the main data.frame and remove
#intermediate ones:
long_ww_lin_w_sample_info <- long_ww_lin_3
rm(long_ww_lin_3, aliases, alias_vector)


###Identify named variants (e.g. VOCs and VUMs)
################################################################################

#Import table of names:
named_variants <- read.csv("../input/voc_voi_vum.csv")

##TO DO: get the variant names to not include "NA" in the pasted string!
#unite(named_variants[2:3], "name_plus_alias", sep = " \\/ ", na.rm = TRUE)

#Default to "other" (unnamed variants will keep this designation):
long_ww_lin_w_sample_info$named_variant_id <- "other"

#Add names for non-recombinant variants:
for (x in seq_len(nrow(long_ww_lin_w_sample_info))) {
  for (y in seq_len(nrow(named_variants))) {
    if (long_ww_lin_w_sample_info$aliases_removed[[x]] ==
        named_variants$lineage_id[[y]] |
        startsWith(long_ww_lin_w_sample_info$aliases_removed[[x]],
                   paste(named_variants$lineage_id[[y]], ".", sep = ""))) {
      message(paste(long_ww_lin_w_sample_info$aliases_removed[x], "->",
                    named_variants$variant_name[y], named_variants$alias_id[y]))
      long_ww_lin_w_sample_info$named_variant_id[x] <-
        paste(named_variants$variant_name[y], named_variants$alias_id[y])
    }
  }
}
rm(named_variants)

#Add names for recombinant variants:
for (x in seq_len(nrow(long_ww_lin_w_sample_info))) {
  lin = long_ww_lin_w_sample_info$aliases_removed[x]
  if (str_starts(lin, "XBB\\.1\\.9\\.2(?![:digit:])")) {
    message(paste(lin, "->", "XBB.1.9.2"))
    long_ww_lin_w_sample_info$named_variant_id[x] <- "XBB.1.9.2"
  } else if (str_starts(lin, "XBB\\.1\\.5(?![:digit:])")) {
    message(paste(lin, "->", "XBB.1.5"))
    long_ww_lin_w_sample_info$named_variant_id[x] <- "XBB.1.5"
  } else if (str_starts(lin, "XBB(?![:alpha:])")) {
    message(paste(lin, "->",  "XBB Other"))
    long_ww_lin_w_sample_info$named_variant_id[x] <- "Other XBB lineage"
  } else if (str_starts(lin, "XBF(?![:alpha:])")) {
    message(paste(lin, "->",  "XBF"))
    long_ww_lin_w_sample_info$named_variant_id[x] <- "XBF"
  } else if (str_starts(lin, "XD(?![:alpha:])")) {
    message(paste(lin, "->",  "XD"))
    long_ww_lin_w_sample_info$named_variant_id[x] <- "XD"
  } else if (str_starts(lin, "XE(?![:alpha:])")) {
    message(paste(lin, "->",  "XE"))
    long_ww_lin_w_sample_info$named_variant_id[x] <- "XE"
  } else if (str_starts(lin, "XF(?![:alpha:])")) {
    message(paste(lin, "->",  "XF"))
    long_ww_lin_w_sample_info$named_variant_id[x] <- "XF"
  } else if (str_starts(lin, "X")) {
    message(paste(lin, "->",  "Other non-XBB recombinant"))
    long_ww_lin_w_sample_info$named_variant_id[x] <- "Other non-XBB recombinant"
  } 
  }

###Export formatted data as csv
################################################################################
write_csv(long_ww_lin_w_sample_info,
          "../cleaned-formatted-data/sequencing-results-formatted-cleaned.csv",
          progress = TRUE, num_threads = 3)
