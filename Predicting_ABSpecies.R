# Load packages
library(caret) # currently have version caret_6.0-62.tar installed (2015-11-23 17:35	3.9M)
library(randomForest)
library(tidyverse)
library(knitr)
library(reshape2)
??melt
########################################################################################

data_in <- read.csv("Sewage_lagoon_all.csv") #call data into R

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
## need to load models from "RFmodel_ABSpecies.RData"
load("RFmodel_ABSpecies.RData")

### run random forest model
dataM <- predict(models4, newdata=data[,2:13], "prob") ##use models to predict bat calls

######################################################################
## Aggregating pulses
M1 <- dataM$rf4

M1$Filename <- filenames_to_keep ##add filenames back in

# this bit will need to be changed depending on how people name their files
M1$Site <- word(M1$Filename, 1, sep=fixed('_'))
M1$Date <- word(M1$Filename, 2, sep=fixed('_'))

######################################################################
# Define functions for computing the category with the highest output probability, and the second and third highest probability. 
# In the case of ties, an NA value is returned:  

get_most_probable <- function(species, rank){
  species <- as.character(species)
  mp <- species[rank == 1];
  if (length(mp) != 1) {mp <- NA};
  mp <- as.character(mp)
  return(mp)
}

get_second_probable <- function(species, rank){
  species <- as.character(species)
  mp <- species[rank == 2];
  if (length(mp) != 1) {mp <- NA};
  mp <- as.character(mp)
  return(mp)
}

get_third_probable <- function(species, rank){
  species <- as.character(species)
  mp <- species[rank == 3];
  if (length(mp) != 1) {mp <- NA};
  mp <- as.character(mp)
  return(mp)
}


bat_data <- M1



