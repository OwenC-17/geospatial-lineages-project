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

remove_zero_qual <- function(map_result, qual) {
  #All other prepocessing must be performed first. Check:
  if (any(nchar(map_result) != nchar(qual))) {
    warning("Number of quality scores does not match number of bases. 
            Preprocessing may be necessary.") 
  }
  
  #The zero-quality character is an exclamation point. Find their locations:
  zero_score_indices <- str_locate_all(qual, "!")[[1]][, 1]
  
  
  if (length(zero_score_indices) > 0) {
    #Split read result into a vector of single characters, then remove those
    #with indices matchin zero-quality scores:
    zeroes_removed <- strsplit(map_result, "")[[1]][-zero_score_indices] %>%
      unlist() %>%
      str_flatten()
    return(zeroes_removed)
  }
  else {
    return(map_result)
  }
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
    mutate(formatted_read_results = mapply(detect_deletions, formatted_read_results)) %>%

    #Remove zero quality bases
    mutate(formatted_read_results = mapply(remove_zero_qual, formatted_read_results, qual_without_gaps)) %>%

    #Remove zeros from quality scores as well:
    mutate(qual_without_gaps = str_remove_all(qual_without_gaps, "\\!")) %>%
  
    #Update coverage so it doesn't include zero-quality bases:
    mutate(coverage_no_zero_qual = str_length(formatted_read_results))
  
  return(prepared_pileup_table)
}


meanPhred <- function(phred_string, offset = 33) {
  mean_quality <- phred_string %>%
    charToRaw() %>%
    as.numeric() %>%
    mean()
  
  return(mean_quality - offset)
}

sdPhred <- function(phred_string) {
  sd_quality <- phred_string %>%
    charToRaw %>%
    as.numeric %>%
    sd()
  
  return(sd_quality)
}

calculate_quality_statistics <- function(formatted_pileup_table) {
  new_table <- formatted_pileup_table %>%
    mutate(mean_quality = mapply(meanPhred, qual_without_gaps)) %>%
    mutate(sd_quality = mapply(sdPhred, qual_without_gaps))
  
  return(new_table)
}




