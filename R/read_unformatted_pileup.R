library(tidyverse)


remove_gaps <- function(formatted_read_results, quality) {
  gap_pos <- str_locate_all(formatted_read_results, "[<>]")[[1]]
  if (length(gap_pos) == 0) {
    qual_string_without_gaps <- quality
  } else {
    gap_pos_vector <- gap_pos[1:(length(gap_pos)/2), 1]
    qual_vector <- strsplit(quality, "")[[1]]
    qual_vector_without_gaps <- qual_vector[-gap_pos_vector]
    qual_string_without_gaps <- paste0(qual_vector_without_gaps, collapse = "")
  }
  return(qual_string_without_gaps)
}

detect_inserts <- function(map_result) {
  temp <- map_result
  #Find 1 or more digits preceded by a plus sign (the digits represent the
  #length of the insert):
  ins_lengths <- unlist(str_extract_all(temp, "(?<=\\+)[:digit:]+"))
  
  #If inserts have been found, ins_lengths will not be null:
  if (!is.null(ins_lengths)) {
    
    #Create a character vector for the bp content of the inserts:
    inserts <- character()
    
    for (ins in ins_lengths) {
      
      #Find the insert sequence starting after the digits:
      ins_start <- str_locate(temp, paste0("\\+", ins))[[2]] + 1
      
      #Use the digits to find the end of the insert sequence:
      ins_end <- ins_start + as.numeric(ins) - 1
      
      #Get the insert sequence as a substring of the whole read_result:
      ins_seq <- str_sub(temp, ins_start, ins_end)
      
      #Remove the insert from the read_result:
      temp <- str_remove(temp, paste0("\\+", ins, ins_seq))
      
      #Add the ins_seq to the inserts vector:
      inserts <- c(inserts, ins_seq)
    }
  }
  return(temp)
}

detect_deletions <- function(map_result) {
  temp <- map_result
  #Find 1 or more digits preceded by a MINUS sign (the digits represent the
  #length of the deletion):
  del_lengths <- unlist(str_extract_all(temp, "(?<=-)[:digit:]+"))
  
  #If deletions have been found, del_lengths will not be null:
  if (!is.null(del_lengths)) {
    
    #Create a character vector for the bp content of the deletions:
    deletions <- character()
    
    for (del in del_lengths) {
      
      #Find the insert sequence starting after the digits:
      del_start <- str_locate(temp, paste0("-", del))[[2]] + 1
      
      #Use the digits to find the end of the insert sequence:
      del_end <- del_start + as.numeric(del) - 1
      
      #Get the insert sequence as a substring of the whole read_result:
      del_seq <- str_sub(temp, del_start, del_end)
      
      #Remove the insert from the read_result:
      temp <- str_remove(temp, paste0("-", del, del_seq))
      
      #Add the ins_seq to the inserts vector:
      deletions <- c(deletions, del_seq)
    }
  }
  return(temp)
}

import_pileup <- function(unformatted_pileup_table) {
  prepared_pileup_table <- unformatted_pileup_table %>%
  
    #Remove extraneous chars ($ = line end, ^. = line start + read score):
    mutate(formatted_read_results = str_remove_all(read_results, "\\$|\\^.")) %>%
    
    #Remove gaps (< or >) and corresponding quality score:
    mutate(qual_without_gaps = mapply(remove_gaps, formatted_read_results, quality)) %>%
    mutate(formatted_read_results = str_remove_all(formatted_read_results, "[<>]")) %>%
    
    #Remove inserts (and store?)
    mutate(formatted_read_results = mapply(detect_inserts, formatted_read_results)) %>%
    
    #Remove deletions
    mutate(formatted_read_results = mapply(detect_deletions, formatted_read_results))
  #Remove zero quality bases

  return(prepared_pileup_table)
}


import_test_1 <- import_pileup(pileup_table)

import_test_2 <- import_pileup(pileup_table2)

import_test_2 <- import_pileup(pileup_test_examples)

import_test_3 <- import_pileup(pileup_table3)

problematic_rows_2 <- import_test_2[str_length(import_test_2$formatted_read_results) != str_length(import_test_2$qual_without_gaps),]
problematic_rows_3 <- import_test_3[str_length(import_test_3$formatted_read_results) != str_length(import_test_3$qual_without_gaps),]
problematic_rows_3b <- import_test_3[str_length(import_test_3$formatted_read_results) != import_test_3$coverage,] 

problematic_rows$diffs <- str_length(problematic_rows$bases_only) - problematic_rows$count

del_example_2 <- as.character("!,,,,,,,,,,-12acgfacfgafgcnnnn.......-3acd..")

dels <- data.frame(str_match_all(del_example_2, "-(\\d+)([:alpha:]+)"))

dels$name <- paste0("DEL_", substr(dels$X3, 1, dels$X2))
dels

table(dels$name)

substr(del_example_2, start = 1, stop = 10)

