sc2_reference_genome <- readLines("../input/MN908947.3.fasta")

sc2_reference_string <- sc2_reference_genome[2:length(sc2_reference_genome)] %>%
  paste0(collapse = "")

nchar(sc2_reference_string)

bases <- paste0(rep(1:nchar(sc2_reference_string), each = 5), c("A", "C", "G", "T", "del"))
head(bases, 25)
tail(bases, 25)

big_table <- matrix(0, nrow = 1, ncol = length(bases), dimnames = list("sample", bases))


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

