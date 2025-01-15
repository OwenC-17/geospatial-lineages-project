library("Rsamtools")
library(tidyverse)

trypileup <- PileupFiles("../input/pileups/BL2_S96.pileup")

fl <- system.file("extdata", "pileup.txt", package="Rsamtools",
                  mustWork=TRUE)

res <- readPileup(fl)

pileup_table <- read_tsv("../input/pileups/BL2_S96.pileup", 
                         col_names = c("ref", "position", "refbase", 
                                       "coverage", "read_results", "quality"))

pileup_table2 <- read_tsv("../input/pileups/105300_S75.pileup", 
                         col_names = c("ref", "position", "refbase", 
                                       "coverage", "read_results", "quality"),
                         col_types = "cicicc")

sum(str_detect(pileup_table2$read_results, "<"))

pileup_table2[str_detect(pileup_table2$read_results, "\\+"),]$read_results

#sum(str_detect(pileup_table2$read_results, "\\.+.*,+"))

original_example_qual <- pileup_table2$quality[148]
raw_example_qual <- charToRaw(pileup_table2$quality[148])
numeric_example_qual <- as.numeric(raw_example_qual)

bases_for_example_qual <- pileup_table2$read_results[148]
str_locate(bases_for_example_qual, "\\.,")

fwrv_point <- substr(original_example_qual, 1815, 1836)
raw_fwrv_point <- charToRaw(fwrv_point)
numeric_fwrv_point <- as.numeric(raw_fwrv_point)
numeric_fwrv_point


pileup_table2$mean_qual <- unlist(test_map_means)

pileup_table2$mean33qual <- pileup_table2$mean_qual - 33
pileup_table2$numeric_qual <- map(pileup_table2$quality, charToRaw) %>%
  map(as.numeric)

pileup_table2$min33qual <- map(pileup_table2$numeric_qual, min) %>%
  unlist() - 33

pileup_table2$max33qual <- map(pileup_table2$numeric_qual, max) %>%
  unlist() - 33

test_map_means <- map(pileup_table2$quality, charToRaw) %>% 
  map(as.numeric) %>%
  map(mean)
head(test_map_means)


lowest_qual_pos <- pileup_table2[pileup_table2$mean_qual == 45, ]


all_qual_characters <- map(pileup_table2$quality, strsplit, split = "") %>% unlist() %>% table()
all_qual_characters

map(all_qual_characters, charToRaw) %>% unlist() %>% as.numeric() - 33

qualstring_a <- pileup_table2$quality[132]
seqstring_a <- pileup_table2$read_results[132]

base_x_qual <- data.frame(base = unlist(strsplit(seqstring_a, "")), qual = unlist(strsplit(qualstring_a, "")))
base_x_qual$rawqual <- map(base_x_qual$qual, charToRaw)
base_x_qual$numqual <- map(base_x_qual$rawqual, as.numeric)




str_extract_all("aaaaaaaaaaaaaaaaaa", "aa")
str_detect("AAAAAAAAAAaaaaaaaa", "aa")


ins_rows <- pileup_table2 %>%
  filter(str_detect(.$read_results, "\\+"))

detect_inserts <- function(map_result) {
  temp <- map_result
  ins_lengths <- unlist(str_extract_all(temp, "(?<=\\+)[:digit:]+"))
  if (!is.null(ins_lengths)) {
    inserts <- character()
    for (ins in ins_lengths) {
      ins_start <- str_locate(temp, paste0("\\+", ins))[[2]] + 1
      ins_end <- ins_start + as.numeric(ins) - 1
      ins_seq <- str_sub(temp, ins_start, ins_end)
      temp <- str_remove(temp, paste0("\\+", ins, ins_seq))
      inserts <- c(inserts, ins_seq)
    }
  }
  return(list(inserts, temp))
}

ins_lengths <- unlist(str_extract_all(test_result_1, "(?<=\\+)[:digit:]+"))

for(ins in ins_lengths) {
  pos <- str_locate(test_result_1, paste0("\\+", ins))[[2]] + 1
  print(pos)
  
  length <- as.numeric(ins)
  
  it <- substr(test_result_1, pos, pos + length - 1)
  print(it)
}

test_example_qual <- pileup_table2$quality[148]
test_result_2 <- pileup_table2$read_results[148]

remove_zero_qual <- function(map_result, qual) {
  if (nchar(map_result) != nchar(qual)) {
   stop("Number of quality scores does not match number of bases. Preprocessing may be necessary.") 
  }
  zero_score_indices <- str_locate_all(qual, "!")[[1]][, 1]

  zeroes_removed <- strsplit(map_result, "")[[1]][-zero_score_indices] %>%
    unlist() %>%
    str_flatten()
    
  return(zeroes_removed)
}

