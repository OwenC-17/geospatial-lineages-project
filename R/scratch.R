chicago_sites <- read.csv("../input/chicago_sites.csv", header = FALSE, col.names = c("number", "sites")) %>%
  mutate(sites = gsub('\\"', '', sites))

long_ww_chicago <- filter(long_ww_lin_w_sample_info, 
                          wwtp_name %in% chicago_sites$sites) %>%
  mutate(wwtp_name = gsub('METROPOLITAN ', 'M', wwtp_name))

#list_of_others <- as.data.frame(table(filter(long_ww_lin_w_sample_info, named_variant_id == "other")$full_lineage_id))

without_date <- filter(long_ww_lin_w_sample_info, is.na(sample_collect_date))
has_date <- filter(long_ww_lin_w_sample_info, !is.na(sample_collect_date))

sum(!unique(without_date$sample_id) %in% unique(has_date$sample_id))
unique(without_date$sample_id)


idph_ww_lin_data


weird_subset <- sample_metadata2[sample_metadata2$sample_id == "102620", ]
rm(weird_subset)
weird_subset2 <- na_remove[na_remove$sample_id == "102620", ]
rm(weird_subset2)
na_remove <- sample_metadata2[!is.na(sample_metadata2$sample_id),]
rm(na_remove)
na_metada <- sample_metadata2[is.na(sample_metadata2$sample_id),]
rm(na_metada)

all_temps <- unique(sample_metadata2$sample_arrival_temp)
weird_temp_ind <- str_which(all_temps, pattern = "^-?[:digit:]+\\.?[:digit:]*$", negate = TRUE)

weird_temps <- all_temps[weird_temp_ind]
weird_temps

#Any string of consecutive dots or commas replaced with a single dot:
weird2 <- str_replace_all(weird_temps, "\\.+|,+", ".")
weird2[1] <- "-12.3." #for testing purposes
weird2

#Any characters that are not numeric or dots or minus signs removed
#AND trailing dots removed:
weird3 <- str_remove_all(weird2, "[^[:digit:]\\.-]|(\\.$)")
weird3
data.frame(weird_temps, as.numeric(weird3))

as.numeric(all_temps[-weird_temp_ind])

sample_metadata_exp <- read_csv("../input/sample_metadata_20250103.csv", guess_max = 10000000)

weird_years <- sample_metadata2[sample_metadata2$year < 2020, ]

NEW_WWTPS <- sample_metadata2[!(sample_metadata2$wwtp_name %in% wwtp_names_dictionary$wwtp_name),]

RECOMBINANT_ONLY <- long_ww_lin_w_sample_info[str_starts(long_ww_lin_w_sample_info$full_lineage_id, "X"),]
RECOMBINANT_SUMMARY <- table(RECOMBINANT_ONLY$full_lineage_id, RECOMBINANT_ONLY$named_variant_id) %>% data.frame() %>%
  filter(Freq > 0) %>%
  mutate(beginning = str_sub(Var1, 1, 7))

BEGINNING_RECOMBINANT <- RECOMBINANT_SUMMARY %>%
  group_by(beginning) %>%
  summarize(noccurences = sum(Freq))


library(microbenchmark)
?microbenchmark

microbenchmark(
  base::strsplit(test_result_2, "", fixed = TRUE),
  base::strsplit(test_result_2, "")
)


pileup_table <- read_tsv("../input/pileups/BL2_S96.pileup", 
                         col_names = c("ref", "position", "refbase", 
                                       "coverage", "read_results", "quality"))

pileup_table2 <- read_tsv("../input/pileups/105300_S75.pileup", 
                          col_names = c("ref", "position", "refbase", 
                                        "coverage", "read_results", "quality"),
                          col_types = "cicicc")

pileup_table3 <- read_tsv("../input/pileups/100284_S68.pileup", 
                          col_names = c("ref", "position", "refbase", 
                                        "coverage", "read_results", "quality"),
                          col_types = "cicicc")

import_test_1 <- import_pileup(pileup_table)

import_test_2 <- import_pileup(pileup_table2)
import_test_2 <- import_test_2 %>%
  mutate(conserved = mapply(is_conserved, formatted_read_results))


#import_test_example <- import_pileup(pileup_test_examples)

import_test_3 <- import_pileup(pileup_table3)
import_test_3 <- import_test_3 %>%
  mutate(conserved = mapply(is_conserved, formatted_read_results))

problematic_rows_2 <- import_test_2[str_length(import_test_2$formatted_read_results) != str_length(import_test_2$qual_without_gaps),]
problematic_rows_3 <- import_test_3[str_length(import_test_3$formatted_read_results) != str_length(import_test_3$qual_without_gaps),]
problematic_rows_3b <- import_test_3[str_length(import_test_3$formatted_read_results) != import_test_3$coverage_no_zero_qual,] 
problematic_rows_3c <- import_test_3 %>%
  mutate(count0qual = str_count(quality, "\\!")) %>%
  filter(coverage != str_length(formatted_read_results) + count0qual)

problematic_rows$diffs <- str_length(problematic_rows$bases_only) - problematic_rows$count

del_example_2 <- as.character("!,,,,,,,,,,-12acgfacfgafgcnnnn.......-3acd..")

dels <- data.frame(str_match_all(del_example_2, "-(\\d+)([:alpha:]+)"))

dels$name <- paste0("DEL_", substr(dels$X3, 1, dels$X2))
dels

table(dels$name)

substr(del_example_2, start = 1, stop = 10)

sum(str_count(import_test_3$formatted_read_results, "\\*"))
