import_pileup <- function(pileup_file) {
  
}

pileup_test_1 <- read_tsv("../input/pileups/BL2_S96.pileup",
                   col_names = c("ref_read", "position", "ref_base", "count",
                                 "results", "quality")) %>%
  mutate(prop_read_start = str_count(results, "\\^") / count) %>%
  mutate(prop_read_end = str_count(results, "\\$") / count) %>%
  mutate(bases_only = str_remove_all(results, pattern = "\\^.|\\$"))
assertthat::are_equal(str_length(pileup_test_1$bases_only), pileup_test_1$count)


pileup_test_2 <- read_tsv("../input/pileups/105300_S75.pileup",
                          col_names = c("ref_read", "position", "ref_base", "count",
                                        "results", "quality")) %>%
  mutate(prop_read_start = str_count(results, "\\^") / count) %>%
  mutate(prop_read_end = str_count(results, "\\$") / count) %>%
  mutate(bases_only = str_remove_all(results, pattern = "\\^.|\\$"))
assertthat::are_equal(str_length(pileup_test_2$bases_only), pileup_test_2$count)


problematic_rows <- pileup_test_2[str_length(pileup_test_2$bases_only) != pileup_test_2$count,]
problematic_rows$diffs <- str_length(problematic_rows$bases_only) - problematic_rows$count

del_example_2 <- as.character("!,,,,,,,,,,-12acgfacfgafgcnnnn.......-3acd..")

dels <- data.frame(str_match_all(del_example_2, "-(\\d+)([:alpha:]+)"))

dels$name <- paste0("DEL_", substr(dels$X3, 1, dels$X2))
dels

table(dels$name)

substr(del_example_2, start = 1, stop = 10)

