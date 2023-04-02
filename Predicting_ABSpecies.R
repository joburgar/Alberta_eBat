# Load packages
library(caret) # currently have version caret_6.0-62.tar installed (2015-11-23 17:35	3.9M)
library(randomForest)
library(tidyverse)
library(knitr)
library(reshape2)

source("compute_bat_counts.R")
load("RFmodel_ABSpecies.RData")

########################################################################################

create_Alberta_eBat_output <- function(input_data=input_data, input_folder="NABat_2022_raw/"){
  
  data_in <- read.csv(paste0(input_folder,input_data)) #call data into R
  
  head(data_in)
  # subset to parameters from AnaLook
  data_sub <- data_in[c("FileName","Dur.ms","Fmax.kHz","Fmin.kHz","Fmean.kHz","Tk.ms","Fk.kHz","Quality..","Tc.ms","Fc.kHz","S1.OPS","Sc.OPS")]
  
  # update colnames to reflect AnaLook names
  colnames(data_sub) <- c("Filename","Dur","Fmax","Fmin","Fmean","Tk","Fk","Qk","Tc","Fc","S1","Sc")
  data_sub$Dc <- data_sub$Tc - data_sub$Tk
  data_sub <- data_sub[c("Filename","Dur","Fmax","Fmin","Fmean","Tk","Fk","Qk","Tc","Fc","Dc","S1","Sc")]
  
  # remove rows that are not complete (blank from filtering in AnabatInsight)
  data <- data_sub[complete.cases(data_sub),]
  
  # keep a record of the filenames in the same order as the data
  filenames_to_keep <- data$Filename
  
  ######################################################################
  ### run random forest model
  dataM <- predict(models4, newdata=data[,2:13], "prob") ##use models to predict bat calls
  
  ######################################################################
  ## Aggregating pulses
  M1 <- dataM$rf4
  
  M1$Filename <- filenames_to_keep ##add filenames back in
  
  # this bit will need to be changed depending on how people name their files
  M1$Site <- word(M1$Filename, 1, sep=fixed('_'))
  M1$Date <- word(M1$Filename, 2, sep=fixed('_'))
  glimpse(M1)
  
  ## Create output that mimics output from Alberta eBat
  # exception is that counts file removes all calls with <3 pulses rather than attributing to unkown
  # are kept in summary file though
  
  eBat_output <- compute_bat_counts_fn(bat_data=M1)
  
  return(eBat_output)
}

########################################################################################

rawdatafiles <- list.files("NABat_2022_raw")
input_data <- rawdatafiles[2]

eBat_output <- create_Alberta_eBat_output(input_data = input_data) 
glimpse(eBat_output)

# need to consider naming convention - does it matter for how it's uploaded into annual report? Might just need the counts and summary bits (i.e., grepl)
write.csv(eBat_output$bat_counts,paste0("NABat_2022_output/",str_sub(input_data,1,end=-5),"_counts.csv"), row.names=F)
write.csv(eBat_output$bat_summary,paste0("NABat_2022_output/",str_sub(input_data,1,end=-5),"_summary.csv"), row.names=F)

rm(data_in, dataM, data_sub, data, M1)
