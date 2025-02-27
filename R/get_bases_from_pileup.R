library(tidyverse)
library(bigmemory)

sc2_reference_genome <- readLines("../input/MN908947.3.fasta")

sc2_reference_string <- sc2_reference_genome[2:length(sc2_reference_genome)] %>%
  paste0(collapse = "")

nchar(sc2_reference_string)

bases <- paste0(rep(1:nchar(sc2_reference_string), each = 5), c("A", "C", "G", "T", "del"))
head(bases, 25)
tail(bases, 25)

big_table <- big.matrix(nrow = 1, ncol = length(bases), type = "integer", dimnames = list("sample", bases))


positions = import_test_3$position
refbases = import_test_3$refbase
readresults = import_test_3$formatted_read_results
for (row in 1:nrow(import_test_3)) {
  pos <- positions[row]
  ref <- refbases[row]
  results <- readresults[row]
  
  big_table["sample", paste0(pos, ref)] <- str_count(results, "[,\\.]")
  
  if (ref != "A") {
    big_table["sample", paste0(pos, "A")] <- str_count(results, "[Aa]")
  }
  
  if (ref != "C") {
    big_table["sample", paste0(pos, "C")] <- str_count(results, "[Cc]")
  }
  
  if (ref != "G") {
    big_table["sample", paste0(pos, "G")] <- str_count(results, "[Gg]")
  }
  
  if (ref != "T") {
    big_table["sample", paste0(pos, "T")] <- str_count(results, "[Tt]")
  }
  
  big_table["sample", paste0(pos, "del")] <- str_count(results, "\\*")
}

base_set <- c("A", "C", "G", "T", "\\*")


# Create a reference vector
ref_bases <- c("A", "C", "G", "T")

# Define a function to check each base
check_base <- function(row, base) {
  if (import_test_3$refbase[row] == base) {
    return(str_count(import_test_3$formatted_read_results[row], "[\\.,]"))
  } else {
    return(str_count(toupper(import_test_3$formatted_read_results[row]), base))
  }
}

import_test_3$A <- mapply(check_base, 1:nrow(import_test_3), base = "A")
import_test_3$C <- mapply(check_base, 1:nrow(import_test_3), base = "C")
import_test_3$G <- mapply(check_base, 1:nrow(import_test_3), base = "G")
import_test_3$T <- mapply(check_base, 1:nrow(import_test_3), base = "T")
import_test_3$del <- mapply(check_base, 1:nrow(import_test_3), base = "\\*")

itest_3_bases <- import_test_3 %>%
  pivot_longer(cols = c("A", "C", "G", "T", "del"), names_to = "detected_base")

ggplot(filter(itest_3_bases, conserved == FALSE)[1:500,], aes(x = as.character(position), y = value, fill = detected_base)) + geom_col(position = "fill")

# Apply the function to each row and each base
result <- apply(import_test_3, 1, function(row_idx) {
  sapply(ref_bases, function(base) check_base(as.numeric(row_idx), base))
})

check_base(15000, "T")
