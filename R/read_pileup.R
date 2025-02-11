####Function to import 1 table
import_muts <- function(pileup_file) {
  pileup <- read_tsv(pileup_file) %>% #read in the file
    replace(is.na(.), 0)
  muts <- filter(pileup, Status != "CONSERVED")
  count_table <- data.frame(row.names = "count")
  for (r in 1:nrow(muts)) {
    example <- muts[r,]
    message(paste("Now processing:", example$Position, "\n"))
    if (example$A > 0) {
      mutationA <- paste(example$Ref_Base, example$Position, "A", sep = "")
      dfA <- data.frame(mutationA)
      dfA$mutationA <- example$A
      names(dfA)[1] <- mutationA
      count_table <- cbind(count_table, dfA)
    }
    if (example$C > 0) {
      mutationC <- paste(example$Ref_Base, example$Position, "C", sep = "")
      dfC <- data.frame(mutationC)
      dfC$mutationC <- example$C
      names(dfC)[1] <- mutationC
      count_table <- cbind(count_table, dfC)
    }
    if (example$`T` > 0) {
      mutationT <- paste(example$Ref_Base, example$Position, "T", sep = "")
      dfT <- data.frame(mutationT)
      dfT$mutationT <- example$`T`
      names(dfT)[1] <- mutationT
      count_table <- cbind(count_table, dfT)
    }
    if (example$G > 0) {
      mutationG <- paste(example$Ref_Base, example$Position, "G", sep = "")
      dfG <- data.frame(mutationG)
      dfG$mutationG <- example$G
      names(dfG)[1] <- mutationG
      count_table <- cbind(count_table, dfG)
    }
  }
  return(count_table)
}



whole_table <- data.frame(row.names = "count", "prop")

#######
rm(pileup_muts)
setwd("../input/pileups/Stephanie_New_Pileups_Folder/230503_pileups/")
for (data in list.files()){ #List files in current directory as "data"
  # Create the initial df if none exists yet
  if (!exists("pileup_muts")){ 
    pileup_muts <- import_muts(data) %>%
      mutate(SampleID = data) 
  }
  
  # if df already exists, then append new data 
  if (exists("pileup_muts")){
    temporary <- import_muts(data) %>% #temporary table is necessary so 
      #that only the current sample is included
      mutate(SampleID = data) #assign the file name as sample ID (it will be split later)
    pileup_muts <- bind_rows(pileup_muts, temporary) #dplyr row binding allows 
    #different numbers of columns
    rm(temporary) #save memory <3 
  }
}

##################################



count1 <- import_muts("../input/pileups/test_example_1000rows.tab")
pileup